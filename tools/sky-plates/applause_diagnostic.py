"""Exploratory APPLAUSE same-image controls; retain the preregistered depth gate."""
import argparse
import csv
import json
import math
import pathlib
import statistics as st
from collections import defaultdict
from astropy.io.votable import parse_single_table
from analysis_core import number, separation_arcsec, in_baseline, differential_curve, annual_contrast, eligible_applause


def csv_rows(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def run(out,recon):
    table=parse_single_table(out/"applause_field_cone10.xml").to_table()
    raw=[{k:None if getattr(row[k],"mask",False) else row[k].item() if hasattr(row[k],"item") else row[k] for k in table.colnames} for row in table]
    original=[r for r in csv_rows(recon/"J0754_applause_cone.csv") if r["natmag"]]
    plates={int(r["plate_id"]) for r in original}
    dates={int(r["plate_id"]):float(r["year_mid"]) for r in csv_rows(recon/"applause_exposure_dates.csv") if r["year_mid"]}
    target_id=667401913089980544
    target_catalogue={f'{r["process_id"]}:{r["solution_num"]}':r for r in original}
    def acceptable(r):
        return (number(r.get("phot_range_flags"))==0 and number(r.get("model_prediction"))>=.9
                and 0<number(r.get("natmag_error"))<=.3 and number(r.get("annular_bin"))<=8
                and math.isfinite(number(r.get("natmag"))) and math.isfinite(number(r.get("cat_natmag"))))
    def converted(r):
        return dict(key=f'{r["process_id"]}:{r["solution_num"]}',year=dates[int(r["plate_id"])],
                    mag=number(r["natmag"])-number(r["cat_natmag"]),plate_id=int(r["plate_id"]),
                    process_id=int(r["process_id"]),good=acceptable(r),
                    below_limit=number(r.get("phot_range_flags"))==2)
    grouped=defaultdict(list)
    position={}
    target=[]
    for r in raw:
        if int(r["plate_id"]) not in plates or r.get("gaiaedr3_id") is None:
            continue
        gid=int(r["gaiaedr3_id"])
        if gid==target_id:
            if math.isfinite(number(r["natmag"])):
                target.append(converted(r))
            continue
        if not acceptable(r) or number(r.get("nn_dist"))<20:
            continue
        reference=target_catalogue.get(f'{r["process_id"]}:{r["solution_num"]}')
        if reference is None or not eligible_applause(r,reference):
            continue
        grouped[str(gid)].append(converted(r))
        position[str(gid)]=separation_arcsec(number(r["ra_icrs"]),number(r["dec_icrs"]),118.69125,16.694722222)
    eligible={k:v for k,v in grouped.items() if len({r["plate_id"] for r in v if in_baseline(r["year"])})>=3}
    ids=sorted(eligible,key=lambda k:position[k])[:12]
    controls={k:eligible[k] for k in ids}
    corrected=differential_curve(target,controls,minimum_controls=5)
    # Different reductions/exposures of one physical plate are not independent votes.
    def collapse(rows):
        groups=defaultdict(list)
        for row in rows:groups[row["plate_id"]].append(row)
        return [dict(year=st.median(r["year"] for r in rr),corrected_mag=st.median(r["corrected_mag"] for r in rr),plate_id=pid,good=all(r["good"] for r in rr),below_limit=any(r["below_limit"] for r in rr)) for pid,rr in groups.items()]
    collapsed=collapse(corrected)
    result=dict(scope="Exploratory diagnostic. It does not override the pre-specified photometric-depth cuts.",
                field_rows=len(raw),selected_controls=ids,
                controls={k:dict(rows=len(controls[k]),distance_arcmin=position[k]/60) for k in ids},
                target_rows_with_controls=len(corrected),target_physical_plates=len(collapsed),
                all_rows_diagnostic=annual_contrast(collapsed,minimum_per_year=1),
                valid_rows_only=annual_contrast([r for r in collapsed if r["good"]],minimum_per_year=1),
                points=corrected,plate_collapsed_points=collapsed,
                verdict="INCONCLUSIVE under pre-specified depth gate. Event-epoch catalog values can only be shown as quality-limited diagnostics.")
    (out/"applause_diagnostic.json").write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ["points","plate_collapsed_points","controls"]},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=pathlib.Path,required=True);parser.add_argument("--recon",type=pathlib.Path,required=True)
    args=parser.parse_args();run(args.out,args.recon)
