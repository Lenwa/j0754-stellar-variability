"""Independent Gaia-forced scan-density pilot; validated support, no extrapolation."""
import argparse
import csv
import json
import pathlib
import warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from plate_pixels import aperture, fit_calibration, predict_magnitude, blank_positions, fit_isotonic

TARGET = 667401913089980544
PLATES = [413, 446, 4191, 4413, 4599]


def values(column):
    return np.asarray(np.ma.filled(column, np.nan), dtype=float)


def run(out, recon, estimator="quadratic", orientation="x"):
    warnings.simplefilter("ignore", FITSFixedWarning)
    g = Table.read(out/"gaia_field.csv", format="ascii.csv")
    catalogue = Table.read(out/"catalogue_union.fits")
    with (recon/"J0754_applause_cone.csv").open() as f:
        target_rows = list(csv.DictReader(f))
    target_index = int(np.flatnonzero(g["source_id"] == TARGET)[0])
    fit = fit_calibration if estimator == "quadratic" else fit_isotonic
    prefix = "" if estimator == "quadratic" else "isotonic_"
    if orientation != "x":
        prefix += orientation+"_"
    ra, dec = values(g["ra"]), values(g["dec"])
    pmra, pmdec = values(g["pmra"]), values(g["pmdec"])
    colour, rp = values(g["bp_rp"]), values(g["phot_rp_mean_mag"])
    rr = ((ra[:,None]-ra)*np.cos(np.radians(dec[:,None])))**2+(dec[:,None]-dec)**2
    np.fill_diagonal(rr, np.inf)
    nearest = np.sqrt(rr.min(axis=1))*3600
    results = []
    fig, axes = plt.subplots(5, 3, figsize=(10, 14), layout="constrained")
    for row_number, plate in enumerate(PLATES):
        with fits.open(out/f"plate_{plate}_{orientation}.fits", checksum=True) as hdus:
            image, h = hdus[0].data.astype(float), hdus[0].header.copy()
        w = WCS(h)
        year = h["YEAR-AVG"]
        meta = next(r for r in target_rows if int(r["plate_id"]) == plate
                    and any(int(s["scan_id"]) == int(r["scan_id"]) for s in
                            json.loads((out/("scan_downloads.json" if orientation == "x" else "scan_downloads_y.json")).read_text()) if s["plate_id"] == plate))
        process = int(meta["process_id"])
        coefficient_source = "target matched catalogue row"
        if meta["color_term"]:
            coefficient = float(meta["color_term"])
            normal = float(meta["cat_natmag"])
        else:
            coefficients = []
            for r in catalogue[catalogue["process_id"] == process]:
                if np.ma.is_masked(r["gaiaedr3_id"]) or np.ma.is_masked(r["cat_natmag"]):
                    continue
                found = np.flatnonzero(g["source_id"] == int(r["gaiaedr3_id"]))
                if len(found) == 1:
                    i = int(found[0])
                    if np.isfinite(rp[i]+colour[i]) and colour[i] > .1:
                        coefficients.append((float(r["cat_natmag"])-rp[i])/colour[i])
            if len(coefficients) < 5:
                paired = next(r for r in target_rows if int(r["plate_id"]) == plate and r["color_term"])
                coefficient, normal = float(paired["color_term"]), float(paired["cat_natmag"])
                coefficient_source = "Assumed paired-scan colour response; this scan has no catalogue calibration"
            else:
                if np.std(coefficients) > .002:
                    raise ValueError("Cannot recover a consistent catalogue colour coefficient")
                coefficient = float(np.median(coefficients))
                normal = float(rp[target_index]+coefficient*colour[target_index])
                coefficient_source = f"Recovered from {len(coefficients)} other matched catalogue rows; target uncalibrated"
        natural = rp + coefficient*colour
        epoch_ra = ra+(year-2016)*pmra/(3.6e6*np.cos(np.radians(dec)))
        epoch_dec = dec+(year-2016)*pmdec/3.6e6
        x, y = w.world_to_pixel_values(epoch_ra, epoch_dec)
        # Register using historical centroids of other stars, never target intensity.
        offsets = []
        for r in catalogue[catalogue["process_id"] == process]:
            if np.ma.is_masked(r["gaiaedr3_id"]) or int(r["gaiaedr3_id"]) == TARGET:
                continue
            found = np.flatnonzero(g["source_id"] == int(r["gaiaedr3_id"]))
            if len(found) != 1 or not float(r["model_prediction"]) >= .9:
                continue
            i = int(found[0])
            if not np.isfinite(x[i]+y[i]):
                continue
            ax, ay = map(float, w.world_to_pixel_values(float(r["ra_icrs"]), float(r["dec_icrs"])))
            if np.hypot(ax-x[i], ay-y[i]) < 10:
                offsets.append([ax-x[i], ay-y[i]])
        offsets = np.asarray(offsets)
        shift = np.median(offsets, axis=0) if len(offsets) else np.array([0., 0.])
        astrometric_scatter = float(1.4826*np.median(np.hypot(*(offsets-shift).T))) if len(offsets) else None
        x, y = x+shift[0], y+shift[1]
        selected = np.ma.filled((np.abs(colour-colour[target_index]) <= .3)
                    & (np.abs(natural-normal) <= 1.5) & (nearest >= 20)
                    & (values(g["ruwe"]) < 1.4) & np.isfinite(x+y)
                    & (g["source_id"] != TARGET)
                    & (x > 25) & (x < image.shape[1]-25)
                    & (y > 25) & (y < image.shape[0]-25), False)
        indices = np.flatnonzero(selected)
        provisional = [aperture(image, x[i], y[i], polarity=1)["signal"] for i in indices]
        polarity = 1 if np.median(provisional) > 0 else -1
        static_x, static_y = w.world_to_pixel_values(ra, dec)
        # Stars without proper motions still mask possible occupied blank-sky areas.
        blank_x = np.where(np.isfinite(x), x, static_x+shift[0])
        blank_y = np.where(np.isfinite(y), y, static_y+shift[1])
        blanks = blank_positions(blank_x, blank_y, image.shape[1], image.shape[0], seed=20261108+plate)
        if len(blanks) < 100:
            raise ValueError(f"Insufficient blank sky on plate {plate}; no noise limit can be inferred")
        entry = dict(plate_id=plate, scan_id=int(meta["scan_id"]), year=year,
                     process_id=process, expected_natural_mag=normal, polarity=polarity,
                     catalogue_colour_coefficient=coefficient, coefficient_source=coefficient_source,
                     position_references=len(offsets), position_shift_pixels=shift.tolist(),
                     position_scatter_pixels=astrometric_scatter,
                     target_pixel=[float(x[target_index]), float(y[target_index])],
                     selected_controls=[str(g["source_id"][i]) for i in indices],
                     blank_apertures=len(blanks), apertures={})
        for radius in [4, 6, 8]:
            controls = [dict(source_id=str(g["source_id"][i]), expected_mag=float(natural[i]),
                             colour=float(colour[i]), x=float(x[i]), y=float(y[i]),
                             **aperture(image, x[i], y[i], radius, polarity)) for i in indices]
            signals = np.array([c["signal"] for c in controls])
            mags = np.array([c["expected_mag"] for c in controls])
            model = fit(signals, mags)
            residuals = []
            for j in range(len(controls)):
                keep = np.arange(len(controls)) != j
                held_model = fit(signals[keep], mags[keep])
                prediction = predict_magnitude(held_model, signals[j])
                controls[j]["held_out_error"] = None if prediction is None else prediction-mags[j]
                controls[j]["held_out_status"] = ("scored" if prediction is not None else
                    "nonpositive_signal" if signals[j] <= 0 else
                    "no_valid_held_out_model" if held_model is None else "outside_held_out_support")
                if prediction is not None:
                    residuals.append(float(prediction-mags[j]))
            blank_signals = np.array([aperture(image, bx, by, radius, polarity)["signal"] for bx, by in blanks])
            blank_center = float(np.median(blank_signals))
            blank_sigma = float(1.4826*np.median(np.abs(blank_signals-blank_center)))
            target = aperture(image, x[target_index], y[target_index], radius, polarity)
            measured = predict_magnitude(model, target["signal"])
            lower_signal = target["signal"]-1.96*blank_sigma
            upper_signal = target["signal"]+1.96*blank_sigma
            lo = predict_magnitude(model, upper_signal)
            hi = predict_magnitude(model, lower_signal)
            quantiles = np.quantile(residuals, [.025, .975]).tolist() if residuals else [None, None]
            rms = float(np.sqrt(np.mean(np.array(residuals)**2))) if residuals else None
            median_error = float(np.median(np.abs(residuals))) if residuals else None
            good = (model is not None and len(residuals) >= 12 and median_error <= .25
                    and len(offsets) >= 5 and astrometric_scatter < 2 and len(blanks) >= 100)
            detected = target["signal"]-blank_center > 5*blank_sigma
            delta = None if measured is None else measured-normal
            corroborates = bool(good and detected and delta is not None
                                and delta > max(.2, 3*rms))
            result = dict(controls=controls, model=model, held_out_errors=residuals,
                          held_out_count=len(residuals), held_out_rms=rms,
                          selected_control_count=len(controls),
                          positive_control_count=int((signals > 0).sum()),
                          unscored_controls=[dict(source_id=c["source_id"],reason=c["held_out_status"])
                                             for c in controls if c["held_out_error"] is None],
                          held_out_median_abs_error=median_error,
                          blank_center=blank_center, blank_sigma=blank_sigma,
                          blank_threshold_exceedances=int((blank_signals > blank_center+5*blank_sigma).sum()),
                          blanks_at_least_as_strong_as_target=int((blank_signals >= target["signal"]).sum()),
                          blank_samples=blank_signals.tolist(),
                          target=target, measured_mag=measured, delta_mag=delta,
                          descriptive_interval=None if lo is None or hi is None or not residuals else
                          [lo-quantiles[1], hi-quantiles[0]],
                          blank_noise_only_bright_bound_mag=predict_magnitude(model, max(target["signal"],0)+5*blank_sigma),
                          calibration_gate=bool(good), target_detection_gate=bool(detected),
                          corroborates_under_screening_gate=corroborates)
            entry["apertures"][str(radius)] = result
        results.append(entry)
        # The same density stretch is used for every stamp within a physical plate.
        primary = entry["apertures"]["6"]
        near = sorted(indices, key=lambda i: abs(natural[i]-normal))[:2]
        for column, i in enumerate([target_index, *near]):
            ix, iy = int(round(x[i])), int(round(y[i]))
            patch = image[iy-35:iy+36, ix-35:ix+36]
            local = aperture(image, x[i], y[i], polarity=polarity)
            display = polarity*(patch-local["background"])
            ax = axes[row_number, column]
            scale = max(1, primary["target"]["pixel_noise"])
            ax.imshow(display, origin="lower", cmap="gray", norm=AsinhNorm(linear_width=scale*2, vmin=-scale*4, vmax=scale*35))
            ax.add_patch(plt.Circle((35+x[i]-ix,35+y[i]-iy),6,fill=False,color="#e89b44",lw=1))
            ax.set_xticks([]); ax.set_yticks([])
            name = "J0754" if column == 0 else "Comparison"
            ax.set_title(f"{year:.1f}, {h['PLATENUM'].strip()} | {name}\nExpected natural mag {natural[i]:.2f}", fontsize=10)
    fig.suptitle("Independent Hamburg plates: forced positions and stars of comparable expected brightness\nSame display scale within each row; circles mark measurement apertures",fontsize=12)
    fig.savefig(out/f"{prefix}applause_star_stamps.png",dpi=160);plt.close(fig)
    physical = [r for r in results if 1930 <= r["year"] < 1940]
    count = sum(r["apertures"]["6"]["corroborates_under_screening_gate"]
                and all(r["apertures"][str(s)]["delta_mag"] is not None
                        and r["apertures"][str(s)]["delta_mag"] > 0 for s in [4,8]) for r in physical)
    output = dict(method=f"{estimator}: Independent catalogue-forced empirical photographic-density calibration; see dated protocols and design amendment.",
                  analysis_status="Primary pixel experiment" if estimator == "quadratic" else "Secondary exploratory calibration after primary-model failure; not preregistered confirmation",
                  results=results, corroborating_physical_event_plates=int(count),
                  verdict="CORROBORATES under specified screening gates" if count >= 2 else "INCONCLUSIVE under specified screening gates",
                  limitations=["Descriptive intervals omit some photographic and reference-catalogue systematics.",
                               "Blank-noise-only bounds are not validated five-sigma stellar brightness limits and are not used for corroboration.",
                               "Held-out RMS covers scored controls only; nonpositive and out-of-support controls have explicit censor reasons and may expose missing calibration coverage.",
                               "Empirical blank-aperture counts are descriptive; overlapping apertures and residual unmasked faint sources preclude a simple independent-trials significance.",
                               "No calibration extrapolation; missing predicted values are not nondetections.",
                               "Same-era colour information is unavailable; natural magnitudes assume the target catalogue colour."])
    (out/f"{prefix}pixel_results.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for r in results:
        a=r["apertures"]["6"]
        print(json.dumps(dict(plate=r["plate_id"],year=r["year"],controls=len(r["selected_controls"]),
                             calibration_stars=None if a["model"] is None else a["model"]["stars"],
                             held_out=a["held_out_count"],rms=a["held_out_rms"],delta=a["delta_mag"],
                             calibration_pass=a["calibration_gate"],detection=a["target_detection_gate"],
                             corroborates=a["corroborates_under_screening_gate"])))
    print(output["verdict"])


if __name__ == "__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True)
    p.add_argument("--recon",type=pathlib.Path,required=True)
    p.add_argument("--estimator",choices=["quadratic","isotonic"],default="quadratic")
    p.add_argument("--orientation",choices=["x","y"],default="x")
    args=p.parse_args();run(args.out,args.recon,args.estimator,args.orientation)
