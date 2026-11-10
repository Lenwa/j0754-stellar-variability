"""Within-camera ASAS-SN replication with UTC-day block resampling."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

from summarize_asassn import camera_filter, accepted


def parse_row(raw):
    row = dict(camera=raw["Camera"], band=raw["Filter"], day=raw["UT Date"][:10],
               year=int(raw["UT Date"][:4]), month=int(raw["UT Date"][5:7]))
    for key, native in [("jd","HJD"),("mag","mag"),("mag_err","mag_err"),
                        ("flux","flux(mJy)"),("flux_err","flux_err"),("fwhm","FWHM")]:
        try:
            row[key] = float(raw[native])
        except ValueError:
            row[key] = float("nan")
    reasons = []
    if raw["mag"].strip().startswith(">"):
        reasons.append("upper_limit")
    if not all(math.isfinite(row[k]) for k in ("jd","mag","mag_err","flux","flux_err")):
        reasons.append("nonfinite")
    if not (0 < row["mag_err"] < 99 and row["mag"] < 90 and row["flux"] > 0 and row["flux_err"] > 0):
        reasons.append("sentinel_or_invalid_error_flux")
    if camera_filter(row["camera"]) != row["band"]:
        reasons.append("unknown_or_inconsistent_filter")
    row["accepted"] = not reasons
    row["rejection_reasons"] = reasons
    return row


def annual(rows):
    return [dict(year=y, observations=len(a), days=len(set(r["day"] for r in a)),
                 median=float(np.median([r["mag"] for r in a])))
            for y in sorted(set(r["year"] for r in rows))
            if (a := [r for r in rows if r["year"] == y])]


def contrast(rows, first, last, draws=3000, seed=20261108):
    if len(set((r["camera"], r["band"]) for r in rows)) > 1:
        raise ValueError("A contrast must stay within one camera/filter")
    samples = [[r for r in rows if r["year"] == y] for y in (first, last)]
    blocks = []
    for sample in samples:
        days = defaultdict(list)
        for row in sample:
            days[row["day"]].append(row["mag"])
        blocks.append([np.asarray(days[k]) for k in sorted(days)])
    result = dict(years=[first,last], observations=[len(s) for s in samples],
                  days=[len(b) for b in blocks], status="INSUFFICIENT")
    if min(result["observations"]) < 10 or min(result["days"]) < 10:
        return result
    medians = [float(np.median([r["mag"] for r in s])) for s in samples]
    result.update(status="MEASURED", medians=medians, delta=medians[1]-medians[0])
    if draws:
        rng = np.random.default_rng(seed)
        simulations = []
        for groups in blocks:
            simulations.append(np.array([np.median(np.concatenate(
                [groups[i] for i in rng.integers(0,len(groups),len(groups))])) for _ in range(draws)]))
        result["ci95"] = np.quantile(simulations[1]-simulations[0], [.025,.975]).tolist()
    return result


def signed(result, direction):
    return result["status"] == "MEASURED" and direction*result["delta"] > 0


def run(out, ztf):
    with (out / "asassn_v1_lightcurve.csv").open(newline="") as f:
        raw = list(csv.DictReader(f))
    rows = [parse_row(r) for r in raw]
    good = [r for r in rows if r["accepted"]]
    cameras = []
    for camera in sorted(set(r["camera"] for r in rows)):
        all_camera = [r for r in rows if r["camera"] == camera]
        data = [r for r in good if r["camera"] == camera]
        item = dict(camera=camera, band=all_camera[0]["band"], raw=len(all_camera),
                    accepted=len(data), annual=annual(data), first_date=min(r["day"] for r in all_camera),
                    last_date=max(r["day"] for r in all_camera), contrasts={})
        for label, first, last in [("fading",2018,2022),("recovery",2022,2025)]:
            primary = contrast(data,first,last)
            common = set(r["month"] for r in data if r["year"] == first) & set(r["month"] for r in data if r["year"] == last)
            variants = dict(error_le_0p1=[r for r in data if r["mag_err"] <= .1],
                            error_and_fwhm=[r for r in data if r["mag_err"] <= .1 and 0 < r["fwhm"] <= 2.5],
                            common_months=[r for r in data if r["month"] in common],
                            january_april=[r for r in data if 1 <= r["month"] <= 4],
                            october_december=[r for r in data if r["month"] >= 10])
            primary["sensitivities"] = {k:contrast(v,first,last,draws=500) for k,v in variants.items()}
            item["contrasts"][label] = primary
        cameras.append(item)

    # V2 controls are coadded/differently reduced. Annual changes are diagnostics,
    # not an invented exact-exposure correction of the V1 target measurements.
    index = pq.read_table(out / "asassn_index.parquet").to_pylist()
    identifiers = {r["asas_sn_id"]:str(r["gaia_id"]) for r in index}
    controls = pq.read_table(out / "asassn_lightcurves.parquet").to_pylist()
    control_rows = []
    for r in controls:
        if r["asas_sn_id"] == 214748736751 or not accepted(r):
            continue
        d = Time(r["jd"],format="jd").datetime
        control_rows.append(dict(r, band=camera_filter(r["camera"]), year=d.year,
                                 month=d.month, day=d.strftime("%Y-%m-%d")))
    placebos = []
    for sid in sorted(set(r["asas_sn_id"] for r in control_rows)):
        for camera in sorted(set(r["camera"] for r in control_rows if r["asas_sn_id"] == sid)):
            data = [r for r in control_rows if r["asas_sn_id"] == sid and r["camera"] == camera]
            placebos.append(dict(gaia_id=identifiers[sid], camera=camera, band=camera_filter(camera),
                                 fading=contrast(data,2018,2022,draws=0),
                                 recovery=contrast(data,2022,2025,draws=0)))

    verdicts = {}
    for label, direction in [("fading",1),("recovery",-1)]:
        supporting, contradictory = [], []
        for item in cameras:
            if item["band"] != "g":
                continue
            result = item["contrasts"][label]
            if result["status"] != "MEASURED":
                continue
            ci = result["ci95"]
            strong = ci[0] > 0 if direction == 1 else ci[1] < 0
            sensitivity_agrees = all(signed(s,direction) for s in result["sensitivities"].values() if s["status"] == "MEASURED")
            control_changes = [r[label]["delta"] for r in placebos
                               if r["camera"] == item["camera"] and r[label]["status"] == "MEASURED"]
            result["control_diagnostic"] = dict(sources=len(control_changes), pipeline="Sky Patrol V2, coadded; not identical-frame V1 photometry")
            if control_changes:
                median = float(np.median(control_changes))
                result["control_diagnostic"].update(median_delta=median,
                    range=[min(control_changes),max(control_changes)], target_minus_median_control=result["delta"]-median)
                sensitivity_agrees &= direction*(result["delta"]-median) > 0
            if strong and sensitivity_agrees:
                supporting.append(item["camera"])
            elif (ci[1] < 0 if direction == 1 else ci[0] > 0) or not sensitivity_agrees:
                contradictory.append(item["camera"])
        verdicts[label] = dict(status="CORROBORATED_DIRECTION" if len(supporting) >= 2 and not contradictory else "INCONCLUSIVE",
                               supporting_cameras=supporting, contradictory_cameras=contradictory)

    corroborated = all(v["status"] == "CORROBORATED_DIRECTION" for v in verdicts.values())
    summary = dict(native_rows=len(rows), accepted_rows=len(good), filters=dict(Counter(r["band"] for r in rows)),
                   accepted_filters=dict(Counter(r["band"] for r in good)),
                   rejection_reasons=dict(Counter(s for r in rows for s in r["rejection_reasons"])),
                   first_date=min(r["day"] for r in rows), last_date=max(r["day"] for r in rows),
                   partial_year_note="2026 observations end8November;2026 is a partial year.",
                   primary_specification="Original protocol.json; CSV-specific implementation/sensitivities specified after preliminary medians in replication_analysis_specification.json.",
                   uncertainty="3000 UTC-date block bootstrap draws per primary contrast; 500 for sensitivities. Sampling intervals do not include calibration systematics or inter-day correlation.",
                   missing_flags="V1 export has no explicit image quality flag. Sentinel/limit/error screens and seeing/error sensitivities are used; no brightness clipping.",
                   cameras=cameras, control_placebos=placebos, verdicts=verdicts,
                   conclusion=("Independent ASAS-SN photometry reproduces the direction of the modern ZTF dimming and recovery." if corroborated else "Directional replication is not established by both sets of gates."),
                   limits="The V2 controls use coadded/different reduction, not identical V1 frames;2025control coverage is insufficient. Cause, period and connection to1930s remain unresolved.")
    (out / "asassn_replication.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")

    fig, axes = plt.subplots(2,1,figsize=(11,8),layout="constrained",height_ratios=[1.5,1])
    colors = {"bi":"tab:blue","bm":"tab:green","bq":"tab:red"}
    for item in cameras:
        if item["camera"] not in colors:
            continue
        a = [r for r in item["annual"] if r["observations"] >= 10 and r["days"] >= 10]
        reference = next(r["median"] for r in a if r["year"] == 2022)
        axes[0].plot([r["year"]+.5 for r in a],[r["median"]-reference for r in a],"o-",
                     color=colors[item["camera"]],label=f"ASAS-SN {item['camera']}, g")
    previous = json.loads((ztf / "modern_differential.json").read_text())["annual"]
    ref = next(r["median_relative_mag"] for r in previous if r["year"] == 2022)
    axes[0].plot([r["year"]+.5 for r in previous],[r["median_relative_mag"]-ref for r in previous],"s--",
                 color="black",label="ZTF g differential")
    axes[0].set(title=("J0754: two surveys reproduce fading and partial recovery" if corroborated else "J0754: camera-by-camera replication diagnostic"),ylabel="Annual median relative to each stream's 2022 value (mag)",xlabel="Calendar year (2026 is partial, through 8 November)")
    axes[0].invert_yaxis()
    axes[0].legend(ncol=2,fontsize=9)
    axes[0].grid(alpha=.2)
    for i,camera in enumerate(colors):
        values = [r["fading"]["delta"] for r in placebos if r["camera"] == camera and r["fading"]["status"] == "MEASURED"]
        axes[1].scatter([i]*len(values),values,s=22,color="grey",alpha=.6,label="Comparison stars (V2)" if i == 0 else None)
        r = next(c for c in cameras if c["camera"] == camera)["contrasts"]["fading"]
        axes[1].errorbar(i,r["delta"],yerr=[[r["delta"]-r["ci95"][0]],[r["ci95"][1]-r["delta"]]],
                         fmt="D",color=colors[camera],capsize=6,label="J0754 (new V1), 95% sampling interval" if i == 0 else None)
    axes[1].axhline(0,color="black",lw=.8)
    axes[1].set(xticks=range(3),xticklabels=list(colors),xlabel="ASAS-SN camera",ylabel="2018 â†’ 2022 change (mag)",title="Target fading exceeds the control-star shifts")
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=.2)
    fig.savefig(out / "asassn_replication.png",dpi=160)
    plt.close(fig)
    print(json.dumps(dict(native_rows=len(rows),accepted_rows=len(good),verdicts=verdicts,
        g_contrasts=[dict(camera=c["camera"],fading=c["contrasts"]["fading"],recovery=c["contrasts"]["recovery"])
                     for c in cameras if c["band"] == "g"]),indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--ztf",type=Path,required=True)
    a = p.parse_args()
    run(a.out,a.ztf)
