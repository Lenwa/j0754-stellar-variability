"""Offline, reproducible J0754 pilot; raw downloads remain outside source."""
import argparse
import csv
import datetime as dt
import json
import math
import pathlib
import statistics as st
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_core import (number, separation_arcsec, in_baseline,
                           quality_dasch, differential_curve, annual_contrast)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def decimal_year(jd):
    date = dt.datetime(2000,1,1,12)+dt.timedelta(days=float(jd)-2451545)
    start, end = dt.datetime(date.year,1,1), dt.datetime(date.year+1,1,1)
    return date.year+(date-start).total_seconds()/(end-start).total_seconds()


def key(row):
    return ":".join(str(row.get(k) or "0") for k in ("series","plate_number","mosaic_number","solution_number"))


def reduce_rows(path, exposure_map, scales, strict=False):
    raw = read_csv(path)
    reduced = []
    for row in raw:
        if not quality_dasch(row, strict=strict):
            continue
        metadata = exposure_map.get(key(row), {})
        scale = scales.get(row["series"])
        reduced.append(dict(key=key(row),year=decimal_year(row["date_jd"]),
                            jd=float(row["date_jd"]),mag=float(row["magcal_magdep"]),
                            err=float(row["magcal_local_rms"]),series=row["series"],
                            plate_number=row["plate_number"],
                            colorterm=number(metadata.get("medianColortermApass")),
                            narrow=scale is not None and scale < 400))
    return reduced, dict(total=len(raw),detections=sum(bool(r["magcal_magdep"]) for r in raw),passed=len(reduced))


def annual(rows, field="corrected_mag"):
    groups = defaultdict(list)
    for row in rows:
        groups[math.floor(row["year"])].append(row[field])
    return [(year+.5,st.median(v),len(v)) for year,v in sorted(groups.items()) if len(v)>=3]


def gap_plot(ax, x, y, *args, **kwargs):
    """Do not imply measurements across long gaps in photographic coverage."""
    xx,yy=[],[]
    for index,(a,b) in enumerate(zip(x,y)):
        if index and a-x[index-1]>2.5:
            xx.append(float("nan"));yy.append(float("nan"))
        xx.append(a);yy.append(b)
    return ax.plot(xx,yy,*args,**kwargs)


def normalize(rows, field):
    baseline = [r[field] for r in rows if in_baseline(r["year"])]
    return st.median(baseline) if baseline else 0.0


def write_rows(path, rows):
    if not rows:
        return
    fields = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader();writer.writerows(rows)


def safe_json(value):
    if isinstance(value,float) and not math.isfinite(value):
        return None
    if isinstance(value,dict):
        return {k:safe_json(v) for k,v in value.items()}
    if isinstance(value,list):
        return [safe_json(v) for v in value]
    return value


def run(out, recon):
    exposures = read_csv(out/"target_exposures.csv")
    exposure_map = {}
    for row in exposures:
        k = ":".join(str(row.get(f) or "0") for f in ("series","platenum","mosnum","solnum"))
        exposure_map[k] = row
    series = json.loads((out/"plate_series.json").read_text(encoding="utf-8-sig"))
    scales = {r["series"]:number(r.get("fitted_plate_scale")) for r in series}
    selected = json.loads((out/"selected_controls.json").read_text())
    target, target_counts = reduce_rows(recon/"J0754_apass_lightcurve.csv",exposure_map,scales)
    strict_target, strict_counts = reduce_rows(recon/"J0754_apass_lightcurve.csv",exposure_map,scales,strict=True)
    controls, control_counts = {}, {}
    for idx, chosen in enumerate(selected):
        name=f"control_{idx+1:02}"
        controls[name], control_counts[name] = reduce_rows(out/(name+".csv"),exposure_map,scales)
        control_counts[name].update(star=chosen["ref_text"],color=chosen["color"],catalog_B=chosen["stdmag"])
    corrected = differential_curve(target,controls)
    strict_controls = {f"control_{idx+1:02}":reduce_rows(out/f"control_{idx+1:02}.csv",exposure_map,scales,strict=True)[0] for idx in range(len(selected))}
    strict_corrected = differential_curve(strict_target,strict_controls)
    variants = {
        "primary": corrected,
        "strict_flags":strict_corrected,
        "narrow_only":[r for r in corrected if r["narrow"]],
        "blue_response":[r for r in corrected if -.3 <= r["colorterm"] <= .7],
        "narrow_blue":[r for r in corrected if r["narrow"] and -.3 <= r["colorterm"] <= .7],
    }
    stats={name:dict(points=len(rows),**annual_contrast(rows)) for name,rows in variants.items()}
    stats["raw_quality_filtered"] = annual_contrast(target,value_key="mag")
    per_series = {s:dict(points=sum(r["series"]==s for r in corrected),**annual_contrast([r for r in corrected if r["series"]==s])) for s in sorted(set(r["series"] for r in corrected))}
    sparse_series = {s:annual_contrast([r for r in corrected if r["series"]==s],minimum_per_year=1) for s in per_series}
    leave_series = {s:annual_contrast([r for r in corrected if r["series"]!=s]) for s in per_series}
    placebo_rows, placebo_stats = {}, {}
    for name, rows in controls.items():
        placebo = differential_curve(rows,{k:v for k,v in controls.items() if k!=name})
        placebo_rows[name] = placebo
        placebo_stats[name] = annual_contrast(placebo)
    # APPLAUSE counts do not imply valid photometry at the event epoch.
    dates = read_csv(recon/"applause_exposure_dates.csv")
    date_map = {r["plate_id"]:float(r["year_mid"]) for r in dates if r["year_mid"]}
    ap = []
    for row in read_csv(recon/"J0754_applause_cone.csv"):
        if not row["natmag"]:
            continue
        year=date_map[row["plate_id"]]
        good = (int(row["phot_range_flags"])==0 and 0<number(row["natmag_error"])<=.3
                and number(row["model_prediction"])>=.9 and int(row["annular_bin"])<=8)
        ap.append(dict(year=year,mag=number(row["natmag"])-number(row["cat_natmag"]),
                       err=number(row["natmag_error"]),plate_id=row["plate_id"],
                       process_id=row["process_id"],colorterm=number(row["color_term"]),
                       below_limit=int(row["phot_range_flags"])==2,good=good))
    ap_event=[r for r in ap if 1930<=r["year"]<1940]
    ap_valid=[r for r in ap if r["good"]]
    ap_summary={"all_calibrated_rows":len(ap),"event_rows":len(ap_event),
                "event_plates":len(set(r["plate_id"] for r in ap_event)),
                "event_below_limit":sum(r["below_limit"] for r in ap_event),
                "valid_event_rows":sum(r["good"] for r in ap_event),
                "valid_rows":len(ap_valid),
                "conclusion":"No event-epoch rows meet the pre-specified depth/quality requirements; pipeline photometry cannot provide clean independent confirmation."}
    catalogue=read_csv(out/"comparison_catalog_20arcmin.csv")
    nearest=[]
    for row in catalogue:
        distance=separation_arcsec(number(row["ra_deg"]),number(row["dec_deg"]),118.69125,16.694722222)
        if distance<60:
            nearest.append(dict(star=row["ref_text"],distance_arcsec=distance,B=row["stdmag"],color=row["color"]))
    # Yearly, rather than pointwise, intervals prevent dense years setting the weight.
    result=dict(target="DASCH J075445.9+164141",protocol="protocol.json",target_quality=target_counts,
                strict_target_quality=strict_counts,comparison_stars=control_counts,
                differential_points=len(corrected),minimum_same_plate_controls=5,
                median_same_plate_controls=st.median([r["controls"] for r in corrected]) if corrected else None,
                variants=stats,series=per_series,leave_one_series_out=leave_series,
                exploratory_sparse_series=sparse_series,
                placebo_controls=placebo_stats,applause=ap_summary,nearby_catalogue_entries=nearest,
                assumptions=["All selected controls retained; robust ensemble medians rather than target-based control pruning.",
                             "Annual medians require at least three measurements; each year receives one vote.",
                             "Bootstrap resamples years independently and does not establish a discovery significance or cover every correlated systematic.",
                             "No brightness sigma-clipping or fitting of a physical mechanism.",
                             "No non-detection was converted into a measured magnitude.",
                             "exploratory_sparse_series relaxes the annual minimum to one point to inspect within-series overlap; it is a secondary, less conservative diagnostic and does not replace the primary calculation."])
    (out/"pilot_results.json").write_text(json.dumps(safe_json(result),indent=2,allow_nan=False),encoding="utf-8")
    write_rows(out/"dasch_quality_filtered.csv",target)
    write_rows(out/"dasch_differential.csv",corrected)
    write_rows(out/"applause_target_quality.csv",ap)
    for name,values in placebo_rows.items():
        write_rows(out/(name+"_placebo.csv"),values)
    make_plot(out,target,corrected,placebo_rows,ap,stats)
    print(json.dumps(dict(target_quality=target_counts,differential_points=len(corrected),variants=stats,applause=ap_summary),indent=2))


def make_plot(out,target,corrected,placebos,ap,stats):
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
    fig, axes=plt.subplots(3,1,figsize=(12,10),sharex=True,layout="constrained")
    ax=axes[0]
    ax.scatter([r["year"] for r in target],[r["mag"] for r in target],s=10,alpha=.2,color="#6686a0",label="Individual observations passing cuts")
    y=annual(target,"mag")
    if y: gap_plot(ax,[v[0] for v in y],[v[1] for v in y],"o-",ms=4,color="#193f5b",label="Year medians (3+ observations; gaps left open)")
    ax.set_ylabel("DASCH B magnitude\nfainter â†“");ax.invert_yaxis();ax.legend(loc="upper left",fontsize=9)
    ax.set_title("J0754: checking the reported 1930s dimming",loc="left",fontweight="bold",fontsize=16)
    ax=axes[1]
    for rows in placebos.values():
        center=normalize(rows,"corrected_mag");y=annual(rows)
        if y: gap_plot(ax,[v[0] for v in y],[v[1]-center for v in y],color="#aeb6bd",alpha=.6,lw=.8)
    center=normalize(corrected,"corrected_mag");y=annual(corrected)
    if y: gap_plot(ax,[v[0] for v in y],[v[1]-center for v in y],"o-",color="#b44732",ms=4,label="Target after same-plate ensemble correction")
    ax.plot([],[],color="#aeb6bd",label="12 controls treated as mock targets")
    ax.axhline(0,color="#555",lw=.6);ax.set_ylabel("Relative B magnitude\nfainter â†“");ax.invert_yaxis();ax.legend(loc="upper left",fontsize=9)
    ax=axes[2]
    good=[r for r in ap if r["good"]];bad=[r for r in ap if not r["good"]]
    for rows,color,label,filled in [(good,"#276782","Pass depth / quality cuts",True),(bad,"#b8753e","Fail cuts; event points below plate limit",False)]:
        if rows:
            ax.errorbar([r["year"] for r in rows],[r["mag"] for r in rows],yerr=[r["err"] for r in rows],fmt="o",ms=5,mfc=color if filled else "white",mec=color,ecolor=color,alpha=.8,label=label,capsize=2)
    ax.axhline(0,color="#555",lw=.6);ax.set_ylabel("APPLAUSE natural mag\nminus Gaia synthetic mag");ax.invert_yaxis();ax.legend(loc="upper left",fontsize=9)
    ax.set_xlabel("Year")
    ax.text(.02,.06,"APPLAUSE event measurements fail the depth cuts. They are not independent confirmation.",transform=ax.transAxes,fontsize=10,color="#8c4e24")
    for ax in axes:
        ax.axvspan(1930,1940,color="#e8bd67",alpha=.2,zorder=-3)
        ax.set_xlim(1888,1994);ax.grid(axis="y",alpha=.15)
    fig.savefig(out/"j0754_pilot.png",dpi=180)
    fig.savefig(out/"j0754_pilot.pdf")
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=pathlib.Path,required=True);parser.add_argument("--recon",type=pathlib.Path,required=True)
    args=parser.parse_args();run(args.out.resolve(),args.recon.resolve())
