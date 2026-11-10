"""Bounded fixed-Gaia-period photometry test; see saved pre-fit protocol.

Run from repository root with python -B
tools/sky-plates/orbit_photometry.py.
No downloads, period search, brightness clipping, or source mutation.
"""
import csv
import json
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from asassn_replication import parse_row

PERIOD = 35.9204566
T0 = 2457000.
OUT = Path("results/sky-plates-j0754-models-20261108/short_period")
CHECKS = Path("results/sky-plates-j0754-checks-20261108")


def daily(rows):
    buckets = defaultdict(list)
    for r in rows:
        buckets[(r["group"], r.get("day", int(r["t"]-.5)))].append(r)
    result = []
    for key in sorted(buckets):
        a = buckets[key]
        result.append(dict(t=float(np.mean([r["t"] for r in a])),
                           y=float(np.median([r["y"] for r in a])),
                           err=float(np.sqrt(sum(r["err"]**2 for r in a))/len(a)),
                           group=key[0], day=str(key[1]), count=len(a)))
    return sorted(result, key=lambda r: (r["t"], r["group"]))


def signal(t, period=PERIOD, harmonics=1):
    phase = 2*np.pi*(np.asarray(t)-T0)/period
    return np.column_stack([f(k*phase) for k in range(1, harmonics+1)
                            for f in (np.cos, np.sin)])


def nuisance(t, groups, mode):
    t = np.asarray(t)
    groups = np.asarray(groups)
    columns = [(groups == g).astype(float) for g in np.unique(groups)]
    if mode.startswith("spline"):
        step = int(mode[6:])
        knots = np.arange(np.floor(t.min()/step)*step,
                          np.ceil(t.max()/step)*step+step/2, step)
        columns.extend(np.interp(t, knots, np.eye(len(knots))[j])
                       for j in range(len(knots)))
    elif mode == "season":
        seasons = np.floor((t-2456900)/365.25)
        for s in np.unique(seasons):
            mask = seasons == s
            columns.extend([mask.astype(float),
                            mask*(t-np.mean(t[mask]))/180.])
    elif mode == "group_linear":
        for g in np.unique(groups):
            mask = groups == g
            columns.append(mask*(t-np.mean(t[mask]))/27.)
    return np.column_stack(columns)


def geometry(t, err, groups, mode, harmonics):
    w = 1/np.asarray(err)
    n = nuisance(t, groups, mode)*w[:, None]
    u, d, _ = np.linalg.svd(n, full_matrices=False)
    rank = int(np.sum(d > max(n.shape)*np.finfo(float).eps*d[0]))
    q = u[:, :rank]
    s = signal(t, PERIOD, harmonics)*w[:, None]
    sr = s-q@(q.T@s)
    op = np.linalg.pinv(sr)*w[None, :]
    retained = np.sum(sr**2, axis=0)/np.sum(s**2, axis=0)
    # Rotation-invariant worst/best signal-energy retention, useful when phase
    # happens to align with a sector's fitted slope.
    ss = s.T@s
    sd, su = np.linalg.eigh(ss)
    invroot = (su/np.sqrt(np.maximum(sd, 1e-30)))@su.T
    eigen_retained = np.linalg.eigvalsh(invroot@(sr.T@sr)@invroot)
    return w, q, sr, op, rank, retained, eigen_retained


def fit(t, y, err, groups, mode="spline365", harmonics=1):
    t, y, err = map(np.asarray, (t, y, err))
    w, q, sr, op, rank, retained, eig = geometry(t, err, groups, mode, harmonics)
    yw = y*w
    yr = yw-q@(q.T@yw)
    beta = op@y
    resid = yr-sr@beta
    dof = max(1, len(t)-rank-len(beta))
    scale = float(resid@resid/dof)
    cov = np.linalg.pinv(sr.T@sr)*scale
    null_sse = float(yr@yr)
    chi = float(resid@resid)
    params = []
    for k in range(harmonics):
        b = beta[2*k:2*k+2]
        amplitude = float(np.linalg.norm(b))
        unit = b/max(amplitude, 1e-30)
        params.append(dict(harmonic=k+1, semi_amplitude_mag=amplitude,
            amplitude_white_sigma=float(np.sqrt(max(0, unit@cov[2*k:2*k+2, 2*k:2*k+2]@unit))),
            cosine_phase_cycles=float(np.arctan2(b[1], b[0])/(2*np.pi)%1)))
    return dict(n=len(t), start_jd=float(min(t)), end_jd=float(max(t)),
        mode=mode, harmonics=harmonics, beta=beta.tolist(), covariance_white=cov.tolist(),
        parameters=params, nuisance_rank=rank, nuisance_columns=nuisance(t, groups, mode).shape[1],
        signal_condition=float(np.linalg.cond(sr)), signal_energy_retained=retained.tolist(),
        phase_energy_retained_range=eig.tolist(),
        weighted_sse=chi, null_weighted_sse=null_sse,
        fraction_sse_improvement=1-chi/null_sse if null_sse else 0,
        delta_bic=float(len(beta)*np.log(len(t))+len(t)*np.log(max(chi,1e-300)/max(null_sse,1e-300))),
        residual_rms_mag=float(np.sqrt(np.mean((resid/w)**2))),
        _residual=(resid/w), _null_residual=yr/w, _op=op,
        _periodic=signal(t, PERIOD, harmonics)@beta,
        _nuisance=(yw-resid-signal(t, PERIOD, harmonics)@beta*w)/w)


def heldout(t, y, err, groups, mode, beta):
    w, q, sr, _, rank, _, _ = geometry(t, err, groups, mode, len(beta)//2)
    yr = np.asarray(y)*w
    yr -= q@(q.T@yr)
    null = float(yr@yr)
    e = yr-sr@np.asarray(beta)
    alternate = float(e@e)
    return dict(n=len(t), fraction_sse_improvement=1-alternate/null if null else 0,
                weighted_sse_null=null, weighted_sse_frozen=alternate,
                white_noise_scaled_delta_chi2=(null-alternate)/(null/max(1,len(t)-rank)),
                frozen_coefficients=list(beta))


def arrays(rows, floor=.02, equal=False):
    return (np.array([r["t"] for r in rows]), np.array([r["y"] for r in rows]),
            np.ones(len(rows)) if equal else np.maximum(floor, [r["err"] for r in rows]),
            np.array([r["group"] for r in rows]))


def clean(item):
    if isinstance(item, dict):
        return {k: clean(v) for k, v in item.items() if not k.startswith("_")}
    if isinstance(item, list):
        return [clean(v) for v in item]
    if isinstance(item, np.ndarray):
        return item.tolist()
    return item


def injections(t, y, err, groups, mode, fitted, draws=200):
    """Injection against real baseline and circular-shift residual sensitivity.

    Circular shifts preserve order per camera, but irregular cadence means this
    is a sensitivity experiment, not calibrated correlated-noise significance.
    """
    op = fitted["_op"]
    s = signal(t)
    baseline = np.asarray(fitted["beta"])
    rng = np.random.default_rng(20261108)
    res = fitted["_null_residual"]
    simulated = []
    for _ in range(draws):
        shifted = np.empty(len(t))
        for g in np.unique(groups):
            mask = groups == g
            shifted[mask] = np.roll(res[mask], rng.integers(0, sum(mask)))
        simulated.append(op@shifted)
    sim = np.array(simulated)
    items = []
    for amp in (.01, .03, .10):
        paired_errors, total_amp, phase_errors, simulation_amp = [], [], [], []
        for phase in np.arange(12)*2*np.pi/12:
            injected = amp*np.array([np.cos(phase), np.sin(phase)])
            recovered = op@(y+s@injected)
            paired_errors.append(float(np.linalg.norm(recovered-baseline-injected)))
            total_amp.append(float(np.linalg.norm(recovered)))
            phase_errors.append(float(np.angle(complex(*recovered)/complex(*injected))/(2*np.pi)))
            simulation_amp.extend(np.linalg.norm(sim+injected, axis=1))
        items.append(dict(injected_semi_amplitude_mag=amp, phases=12,
            maximum_paired_coefficient_error_mag=max(paired_errors),
            recovered_total_amplitude_range_mag=[min(total_amp),max(total_amp)],
            maximum_total_phase_error_cycles=max(abs(np.array(phase_errors))),
            circular_shift_recovery_amplitude_p05_p50_p95=np.quantile(simulation_amp,[.05,.5,.95]).tolist()))
    return dict(paired_recovery_note="Joint fitting has unit linear response where identifiable; this does not prove the baseline contains a stable signal.",
        circular_shift_draws=draws, circular_shift_null_amplitude_p50_p95=np.quantile(np.linalg.norm(sim,axis=1),[.5,.95]).tolist(),
        circular_shift_limit="Within-camera residual circular shifts preserve sample order but distort irregular time spacings; diagnostic sensitivity, not calibrated false alarm probability.",
        cases=items)


def epoch_test(rows, mode, harmonics=1, equal=False, boundary=2459580.5):
    a = arrays(rows, equal=equal)
    full = fit(*a, mode, harmonics)
    train_mask = a[0] < boundary
    if min(sum(train_mask),sum(~train_mask)) < 12:
        return dict(full=full, status="INSUFFICIENT_INDEPENDENT_EPOCHS")
    train_args = tuple(x[train_mask] for x in a)
    val_args = tuple(x[~train_mask] for x in a)
    train = fit(*train_args, mode, harmonics)
    val = fit(*val_args, mode, harmonics)
    score = heldout(*val_args, mode, train["beta"])
    b = np.array(train["beta"])-val["beta"]
    c = np.array(train["covariance_white"])+val["covariance_white"]
    phase_difference = ((val["parameters"][0]["cosine_phase_cycles"]-
                        train["parameters"][0]["cosine_phase_cycles"]+.5)%1)-.5
    return dict(full=full, train=train, validation=val, heldout=score,
        phase_difference_cycles=phase_difference,
        coefficient_discrepancy_white_chi2=float(b@np.linalg.pinv(c)@b),
        status="DIAGNOSTIC_ONLY_WHITE_UNCERTAINTY")


def load_asas():
    path = Path("results/sky-plates-j0754-independent-20261108/asassn_v1_lightcurve.csv")
    with path.open(newline="") as f:
        parsed = [parse_row(r) for r in csv.DictReader(f)]
    grouped = {}
    for band in ("g", "V"):
        grouped[band] = daily([dict(t=r["jd"],y=r["mag"],err=r["mag_err"],
            group=r["camera"],day=r["day"]) for r in parsed
            if r["accepted"] and r["band"] == band])
    return grouped, dict(raw=len(parsed), valid=sum(r["accepted"] for r in parsed),
        rejects=dict(Counter(x for r in parsed for x in r["rejection_reasons"])))


def load_tess():
    series = defaultdict(list)
    inventory = []
    for path in sorted(CHECKS.glob("*tess*fits")):
        with fits.open(path) as h:
            data, header = h[1].data, h[0].header
            sector = str(header.get("SECTOR"))
            provider = "SPOC" if "tess-spoc" in path.name else "QLP"
            for col in ("SAP_FLUX", "PDCSAP_FLUX", "KSPSAP_FLUX", "DET_FLUX"):
                if col not in data.names:
                    continue
                key = provider+"_"+col
                flux = np.asarray(data[col], dtype=float)
                t = np.asarray(data["TIME"], dtype=float)+2457000
                mask = (np.asarray(data["QUALITY"]) == 0)&np.isfinite(t)&np.isfinite(flux)&(flux>0)
                errcol = col+"_ERR"
                error = np.asarray(data[errcol], dtype=float) if errcol in data.names else None
                if error is not None:
                    mask &= np.isfinite(error)&(error>0)
                ref = np.median(flux[mask])
                for i in np.flatnonzero(mask):
                    series[key].append(dict(t=float(t[i]), y=float(-2.5*np.log10(flux[i]/ref)),
                        err=float(2.5/np.log(10)*error[i]/flux[i]) if error is not None else .001,
                        group=sector))
                inventory.append(dict(file=path.name, sector=sector, series=key,
                    quality_zero_points=int(sum(mask)), median_flux=float(ref),
                    start_jd=float(min(t[mask])),end_jd=float(max(t[mask]))))
    return {k: daily(v) for k,v in series.items()}, inventory


def main():
    if not (OUT/"protocol.json").exists():
        raise RuntimeError("Save the prediction/protocol before fitting")
    asas, inv = load_asas()
    g = asas["g"]
    variants = {}
    for mode in ("spline365", "spline180", "spline730", "season"):
        variants[mode] = epoch_test(g, mode)
    variants["equal_weight_spline365"] = epoch_test(g, "spline365", equal=True)
    variants["spline365_plus_harmonic"] = epoch_test(g, "spline365", 2)
    primary = variants["spline365"]
    injection = injections(*arrays(g), "spline365", primary["full"])
    vb = epoch_test(asas["V"], "spline365", boundary=2457754.5)
    rawztf = json.loads((CHECKS/"modern_differential.json").read_text())["measurements"]
    z = daily([dict(t=r["mjd"]+2400000.5,y=r["corrected_relative_mag"],
        err=.02,group="ZTF_differential_g") for r in rawztf])
    ztf = {m: epoch_test(z,m) for m in ("spline365","season")}
    tess_rows, inventory = load_tess()
    tess = {}
    for key, rows in tess_rows.items():
        result = epoch_test(rows, "group_linear", boundary=2460000.)
        result["daily_groups"] = len(rows)
        result["sectors"] = sorted(set(r["group"] for r in rows))
        result["by_sector"] = {}
        for sector in result["sectors"]:
            subset = [r for r in rows if r["group"] == sector]
            if len(subset) >= 6:
                result["by_sector"][sector] = fit(*arrays(subset, floor=.001), "group_linear")
        tess[key] = result
    result = dict(period_days=PERIOD, phase_zero_jd=T0, protocol="protocol.json",
        asassn_inventory=inv, daily_counts={k:len(v) for k,v in asas.items()},
        asassn_g=variants, injection=injection, asassn_V=vb, ztf=ztf,
        tess_inventory=inventory,tess=tess,
        uncertainty="Covariances and coefficient discrepancy assume independent daily noise after nuisance fit and estimate residual variance. They exclude red noise, calibration, period/ephemeris uncertainty and detrending uncertainty. Delta BIC is a descriptive independent-Gaussian statistic only.",
        verdict="NOT_ESTABLISHED: the exact35.9204566-day training waveform worsens independent-epoch fits under every tested ASAS-SN nuisance variant; TESS is reduction-dependent and baseline-limited.")
    (OUT/"fixed_period_results.json").write_text(json.dumps(clean(result),indent=2)+"\n")
    with (OUT/"asassn_g_daily.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(g[0]))
        writer.writeheader();writer.writerows(g)
    plot(g, variants, tess_rows, tess)
    report(result)
    print(json.dumps(dict(asas_primary=clean(primary),injection=injection,
                         ztf={k:clean(v) for k,v in ztf.items()}),indent=2))


def plot(g, variants, tess_rows, tess):
    fig, axes = plt.subplots(2,2,figsize=(13,9),layout="constrained")
    t,y,e,groups = arrays(g)
    f=variants["spline365"]["full"]
    axes[0,0].scatter((t-2457000)/365.25+2014.94,y-f["_nuisance"],s=6,alpha=.35)
    axes[0,0].set(xlabel="Approximate year",ylabel="ASAS-SN g minus camera + slow nuisance (mag)",title="Joint 365-day trend and fixed 35.9204566-day signal")
    phase=(t-T0)/PERIOD%1
    for label,mask,color in (("Train: 2017â€“2021",t<2459580.5,"tab:blue"),
                             ("Validation: 2022â€“2026",t>=2459580.5,"tab:orange")):
        subset=[r for r,m in zip(g,mask) if m]
        ff=fit(*arrays(subset),"spline365")
        tt,yy,_,_=arrays(subset)
        axes[0,1].scatter((tt-T0)/PERIOD%1,yy-ff["_nuisance"],s=7,alpha=.25,color=color)
        p=np.linspace(0,1,200)
        b=ff["beta"]
        axes[0,1].plot(p,b[0]*np.cos(2*np.pi*p)+b[1]*np.sin(2*np.pi*p),color=color,label=label)
    axes[0,1].set(xlabel="Fixed Gaia-candidate phase (zero JD 2457000)",ylabel="ASAS-SN g residual (mag)",title="Training waveform does not persist in validation")
    axes[0,1].legend(fontsize=8)
    labels=[]; amps=[];scores=[]
    for k,v in variants.items():
        labels.append(k);amps.append(v["full"]["parameters"][0]["semi_amplitude_mag"])
        scores.append(v.get("heldout",{}).get("fraction_sse_improvement",np.nan)*100)
    x=np.arange(len(labels))
    axes[1,0].bar(x,amps,color="slategrey")
    axes[1,0].set(xticks=x,xticklabels=labels,ylabel="Full-series semi-amplitude (mag)",title="Sensitivity to explicitly joint slow-trend model")
    axes[1,0].tick_params(axis="x",rotation=40,labelsize=8)
    twin=axes[1,0].twinx();twin.plot(x,scores,"ro");twin.axhline(0,color="red",lw=.5)
    twin.set_ylabel("Held-out SSE improvement (%)",color="red")
    for key in ("QLP_SAP_FLUX","QLP_KSPSAP_FLUX","QLP_DET_FLUX","SPOC_SAP_FLUX","SPOC_PDCSAP_FLUX"):
        if key not in tess_rows: continue
        rr=tess_rows[key]
        xx=np.array([r["t"] for r in rr]); yy=np.array([r["y"] for r in rr])
        axes[1,1].scatter((xx-T0)/PERIOD%1,yy,s=9,alpha=.6,label=key)
    axes[1,1].set(xlabel="Fixed Gaia-candidate phase",ylabel="Per-sector relative TESS magnitude",title="TESS reductions kept separate; sector offsets removed")
    axes[1,1].legend(fontsize=7)
    for ax in axes.flat: ax.grid(alpha=.2)
    fig.suptitle("J0754 bounded fixed-period test â€” detrending and daily correlated noise limit inference")
    fig.savefig(OUT/"fixed_period_diagnostic.png",dpi=170)
    plt.close(fig)


def report(result):
    p=result["asassn_g"]["spline365"]
    lines=["# J0754: fixed 35.9204566-day photometric test", "",
      "The tested stable waveform does not persist across independent epochs. Applying the training waveform worsens held-out ASAS-SN and ZTF fits; TESS cannot establish a coherent counterpart from these reductions and sector baselines.", "",
      "The prediction/protocol was saved before fitting. No period scan or brightness clipping was performed. No external data were downloaded.", "",
      "## ASAS-SN g: primary bounded test", "",
      f"{result['asassn_inventory']['valid']}/{result['asassn_inventory']['raw']} valid individual measurements across V+g; {result['daily_counts']['g']} daily camera groups in g.", "",
      "|Joint nuisance|Full semi-amplitude (mag)|Train amplitude|Validation amplitude|Phase difference (cycles)|Held-out SSE improvement|", "|---|---:|---:|---:|---:|---:|"]
    for k,v in result["asassn_g"].items():
        fmt=lambda f:f["parameters"][0]["semi_amplitude_mag"]
        lines.append(f"|{k}|{fmt(v['full']):.5f}|{fmt(v['train']):.5f}|{fmt(v['validation']):.5f}|{v['phase_difference_cycles']:+.3f}|{v['heldout']['fraction_sse_improvement']*100:+.2f}%|")
    lines.extend(["", "Training ends 31 December 2021; validation begins 1 January 2022. Phase is the cosine-maximum convention at JD 2457000. Each validation fit re-estimates nuisance terms; held-out coefficients remain frozen to training. Positive held-out improvement alone is insufficient for a detection. Validation amplitude is consistent with zero, so its point-estimate phase is poorly defined.","",
       f"Primary full-series amplitude {p['full']['parameters'][0]['semi_amplitude_mag']:.5f} Â± {p['full']['parameters'][0]['amplitude_white_sigma']:.5f} mag (white-noise-only 1Ïƒ); validation {p['validation']['parameters'][0]['semi_amplitude_mag']:.5f} Â± {p['validation']['parameters'][0]['amplitude_white_sigma']:.5f} mag. Train/validation coefficient-discrepancy Ï‡Â²={p['coefficient_discrepancy_white_chi2']:.2f} for two coefficients; this is not calibrated for red noise.","",
       "## Actual-cadence injection/sensitivity", "",
       f"The primary nuisance model retains {p['full']['phase_energy_retained_range'][0]:.3f}â€“{p['full']['phase_energy_retained_range'][-1]:.3f} of signal energy over all sine phases. Correct joint fitting recovers paired coefficient increments to numerical precision, while arbitrary detrending before fitting would attenuate the signal.","",
       f"200 within-camera circular shifts of the real null residual yield median/95th-percentile apparent semi-amplitude {result['injection']['circular_shift_null_amplitude_p50_p95']}. These shifts preserve camera sample order but distort irregular time spacings; they are a sensitivity diagnostic, not a false-alarm calculation.","",
       "|Injected semi-amplitude|Recovered total amplitude range over 12 phases|Maximum total phase error (cycles)|", "|---|---|---:|"])
    for v in result["injection"]["cases"]:
        ar=v['recovered_total_amplitude_range_mag']
        lines.append(f"|{v['injected_semi_amplitude_mag']:.3f} mag|{ar[0]:.5f}â€“{ar[1]:.5f} mag|{v['maximum_total_phase_error_cycles']:.3f}|")
    lines.extend(["", "## Independent ZTF g diagnostic", ""])
    for k,v in result["ztf"].items():
        lines.append(f"- {k}: full semi-amplitude {v['full']['parameters'][0]['semi_amplitude_mag']:.5f} mag; held-out SSE improvement {v['heldout']['fraction_sse_improvement']*100:+.2f}%; phase difference {v['phase_difference_cycles']:+.3f} cycles.")
    lines.extend(["", "ZTF uses existing same-exposure differential estimates and equal daily weights; no per-point uncertainty was present in the retained JSON. V-band is a separate diagnostic in structured results, split at 2017 because most V coverage precedes 2022. ASAS-SN V1 has no explicit image-quality flags; the retained parser rejects limits, invalid errors/flux, sentinels and inconsistent camera filters. There is no brightness clipping.","", "## TESS: short-sector ambiguity", "",
      "QUALITY==0, finite positive flux. All eight files were inspected. QLP SAP and corrected KSPSAP/DET flux, SPOC SAP and PDCSAP flux remain separate. Some are duplicate reductions of the same 2023 epochs. Per-sector medians and fitted offsets/slopes remove longterm information. The fixed 36-day signal has less than one cycle per usual sector; short-sector inferred amplitudes can be strongly degenerate with trends. Joint scores use effectively equal daily weights (0.02 mag floor); individual sector diagnostics use a 0.001 mag floor. Sector 1751 has only five valid daily groups from 18â€“22 January 2026 and is insufficient for its own periodic fit.", "",
      "|Reduction|Sectors|Full semi-amplitude|Phase-energy retained range|Independent 2021 vs 2023+ comparison|", "|---|---|---:|---|---|"])
    for k,v in result["tess"].items():
        f=v["full"]
        outcome=(f"held-out improvement {v['heldout']['fraction_sse_improvement']*100:+.1f}%, phase shift {v['phase_difference_cycles']:+.3f}" if "heldout" in v else "insufficient independent early epoch")
        lines.append(f"|{k}|{','.join(v['sectors'])}|{f['parameters'][0]['semi_amplitude_mag']:.5f} mag|{f['phase_energy_retained_range'][0]:.3f}â€“{f['phase_energy_retained_range'][-1]:.3f}|{outcome}|")
    lines.extend(["", "## Interpretation limits", "",
      "A stable signal at the exact Gaia candidate period is not established when epoch phase and held-out benefit disagree or amplitude changes materially with nuisance/reduction choice. All quoted covariance uses independent daily residuals and excludes correlated/systematic noise, calibration drift, and Gaia period/ephemeris uncertainty. The fixed-period failure would not exclude a nearby drifting signal, rotation, an eclipsing system with non-sinusoidal or intermittent variability, or a companion without a photometric signature. Those hypotheses were not searched here.", "",
      "The candidate Gaia orbit fails recommended gold-sample significance/efficiency cuts; even a photometric match would not establish a companion or the cause of the historical dimming. Maintained source is `tools/sky-plates/orbit_photometry.py`; results are local and uncommitted. No online publication."])
    (OUT/"report.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
