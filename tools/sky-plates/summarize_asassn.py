"""Audit ASAS-SN coverage without interpreting absent rows as constant light."""
import argparse
from collections import Counter
import csv
from html.parser import HTMLParser
import json
import math
from pathlib import Path

from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

TARGET = 214748736751
VARIABLE_URL = "https://asas-sn.osu.edu/variables/05a882b3-3477-5b94-b81b-17a186f6b4c1"
V_CAMERAS = set("ba bb bc bd be bf bg bh".split())
G_CAMERAS = set("bA bB cB bC bD bE bF bG bH bi bj bk bl bm bn bo bp bq br bs bt".split())


def camera_filter(camera):
    return "V" if camera in V_CAMERAS else "g" if camera in G_CAMERAS else None


def accepted(row, require_quality=True):
    keys = ("jd", "mag", "mag_err", "flux", "flux_err")
    return (all(math.isfinite(row[k]) for k in keys)
            and 0 < row["mag_err"] < 99 and row["flux"] > 0 and row["flux_err"] > 0
            and (not require_quality or row["quality"] == "G")
            and camera_filter(row["camera"]) is not None)


def annual(rows):
    if not rows:
        return []
    years = np.array([d.year for d in Time([r["jd"] for r in rows], format="jd").datetime])
    mags = np.array([r["mag"] for r in rows])
    result = []
    for year in sorted(set(years)):
        values = mags[years == year]
        result.append(dict(year=int(year), n=len(values), median=float(np.median(values)),
                           adequate_count=len(values) >= 10))
    return result


def coverage(index, rows, target_id=TARGET):
    by_id = {}
    for row in rows:
        by_id.setdefault(row["asas_sn_id"], []).append(row)
    sources = []
    for item in index:
        data = by_id.get(item["asas_sn_id"], [])
        cameras = []
        for camera in sorted(set(r["camera"] for r in data)):
            raw = [r for r in data if r["camera"] == camera]
            good = [r for r in raw if accepted(r)]
            cameras.append(dict(camera=camera, phot_filter=camera_filter(camera), raw=len(raw),
                                good=len(good), quality_counts=dict(Counter(r["quality"] for r in raw)),
                                first_jd=min(r["jd"] for r in raw), last_jd=max(r["jd"] for r in raw),
                                annual=annual(good)))
        sources.append(dict(gaia_id=str(item["gaia_id"]), asas_sn_id=str(item["asas_sn_id"]),
                            raw_rows=len(data), cameras=cameras))
    n = len(by_id.get(target_id, []))
    return dict(matched_catalogue_sources=len(index), sources_with_measurements=len(by_id),
                target_rows=n, target_status="NO_RETURNED_MEASUREMENTS" if n == 0 else "ROWS_RETURNED",
                sources=sources)


class VisibleText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.skip = False
        self.values = []

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.skip = True

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.skip = False

    def handle_data(self, value):
        if not self.skip and value.strip():
            self.values.append(value.strip())


def run(out, ztf):
    index = pq.read_table(out / "asassn_index.parquet").to_pylist()
    rows = pq.read_table(out / "asassn_lightcurves.parquet").to_pylist()
    audit = coverage(index, rows)
    if audit["target_rows"]:
        raise ValueError("New target measurements require the planned per-camera replication analysis; do not reuse the missing-coverage verdict")
    independent_attempts = []
    for name in ("asassn_target_lightcurve.parquet", "asassn_master_cone_lightcurves.parquet"):
        table = pq.read_table(out / name)
        independent_attempts.append(dict(file=name, rows=table.num_rows,
                                        target_rows=sum(x == TARGET for x in table["asas_sn_id"].to_pylist())))
    if any(r["target_rows"] for r in independent_attempts):
        raise ValueError("Supplementary queries now return target data; the missing-coverage summary no longer applies")
    published = json.loads((out / "asassn_variable_target.json").read_text())
    old = [dict(r, jd=r["hjd"]) for r in published["results"]]
    if published["count"] != len(old):
        raise ValueError("The published JSON is paginated or incomplete")
    parser = VisibleText()
    parser.feed((out / "asassn_variable_target.html").read_text())
    values = parser.values
    if "667401913089980544" not in values:
        raise ValueError("Variable-page Gaia identity check failed")
    claims = {k: values[values.index(k)+1] for k in
              ("Mean VMag", "Amplitude", "Period", "Type", "Classification Probability", "Reference")}
    old_good = [r for r in old if accepted(r, require_quality=False)]
    with (out / "asassn_variable_target.csv").open(newline="") as f:
        csv_export_rows = len(list(csv.DictReader(f)))
    old_cameras = []
    for camera in sorted(set(r["camera"] for r in old)):
        raw = [r for r in old if r["camera"] == camera]
        good = [r for r in old_good if r["camera"] == camera]
        old_cameras.append(dict(camera=camera, phot_filter=camera_filter(camera), raw=len(raw),
                                good=len(good), annual=annual(good),
                                first_date=Time(min(r["jd"] for r in raw), format="jd").iso,
                                last_date=Time(max(r["jd"] for r in raw), format="jd").iso,
                                median=float(np.median([r["mag"] for r in good])),
                                robust_scatter=float(1.4826*np.median(np.abs(
                                    np.array([r["mag"] for r in good])-np.median([r["mag"] for r in good]))))))
    summary = dict(question="Direct Sky Patrol V2 coverage audit for J0754 (separate from completed V1 replication)",
                   verdict="INCONCLUSIVE_INSUFFICIENT_TARGET_COVERAGE", bulk_audit=audit,
                   independent_query_attempts=independent_attempts,
                   older_variable_catalogue=dict(url=VARIABLE_URL, published_rows=len(old),
                       csv_export_rows=csv_export_rows,
                       export_note="Native CSV export is header-only in this snapshot; JSON supplies the 222 rows, also embedded on the HTML page.",
                       screening="Finite samples; 0<mag_err<99, positive flux and flux_err. No quality flags supplied.",
                       cameras=old_cameras, catalogue_claims=claims,
                       interpretation="ROT: and 912-day period are catalogue claims, not a mechanism or a period verified here."),
                   interpretation="Zero returned rows are not a flat light curve. Older V data cannot test the 2022 minimum or 2025 recovery.",
                   next_step=("Completed V1 photometry is available. See asassn_replication.json for the independent modern test; this summary only audits direct V2 and older published data." if (out / "asassn_v1_lightcurve.csv").exists() else "Run or resume the V1 coordinate-photometry job recorded in asassn_v1_job.json."))
    (out / "asassn_summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    (out / "asassn_variable_visible_text.txt").write_text("\n".join(values)+"\n", encoding="utf-8")

    # Separate panels and labels: no absolute cross-filter/survey zero-point alignment.
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.6), sharex=True, layout="constrained")
    for camera in sorted(set(r["camera"] for r in old_good)):
        data = [r for r in old_good if r["camera"] == camera]
        x = Time([r["jd"] for r in data], format="jd").decimalyear
        axes[0].scatter(x, [r["mag"] for r in data], s=9, alpha=.5, label=f"camera {camera}, {camera_filter(camera)}")
        ann = [r for r in annual(data) if r["adequate_count"]]
        axes[0].plot([r["year"]+.5 for r in ann], [r["median"] for r in ann], "o-", color="black", label="Annual median (â‰¥10 points)")
    axes[0].set(ylabel="ASAS-SN published V magnitude", title="J0754: older ASAS-SN data stop in 2018")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[0].axvspan(2019, 2026, alpha=.08, color="grey")
    axes[0].text(2022, 12.22, "No target measurements returned\nfor this interval by Sky Patrol V2", ha="center", fontsize=9)
    previous = json.loads((ztf / "modern_differential.json").read_text())
    # This existing result has a single principal stream; fail rather than guess its schema.
    yr = previous["annual"]
    axes[1].plot([r["year"]+.5 for r in yr], [r["median_relative_mag"] for r in yr], "o-", color="tab:orange")
    axes[1].set(ylabel="ZTF g differential magnitude", xlabel="Calendar year",
                title="Earlier ZTF result â€” the direct V2 catalogue download supplies no target counterpart")
    for ax in axes:
        ax.invert_yaxis()
        ax.grid(alpha=.2)
        ax.set_xlim(2012.7, 2026)
    fig.savefig(out / "asassn_coverage.png", dpi=160)
    plt.close(fig)
    print(json.dumps(dict(verdict=summary["verdict"], target_rows=audit["target_rows"],
                         comparison_sources=audit["sources_with_measurements"], older_rows=len(old),
                         catalogue_claims=claims), indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ztf", type=Path, required=True)
    a = p.parse_args()
    run(a.out, a.ztf)
