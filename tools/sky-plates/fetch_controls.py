"""Download outcome-blind comparison stars and independent APPLAUSE plate data."""
import argparse
import csv
import io
import json
import pathlib
import time
import requests
from analysis_core import number, separation_arcsec, eligible_catalog

BASE = "https://api.starglass.cfa.harvard.edu/public/"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save_csv(path, lines):
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")


def post_csv(api, payload):
    response = requests.post(BASE+api, json=payload, timeout=60)
    response.raise_for_status()
    return response.json()


def controls(out, recon):
    out.mkdir(parents=True,exist_ok=True)
    if not (out/"comparison_catalog_20arcmin.csv").exists():
        save_csv(out/"comparison_catalog_20arcmin.csv",post_csv("dasch/dr7/querycat",dict(refcat="apass",ra_deg=118.69125,dec_deg=16.694722222,radius_arcsec=1200)))
    if not (out/"target_exposures.csv").exists():
        save_csv(out/"target_exposures.csv",post_csv("dasch/dr7/queryexps",dict(ra_deg=118.69125,dec_deg=16.694722222)))
    if not (out/"plate_series.json").exists():
        response=requests.get(BASE+"series",timeout=60)
        response.raise_for_status()
        (out/"plate_series.json").write_text(json.dumps(response.json(),indent=2),encoding="utf-8")
    catalog = rows(out/"comparison_catalog_20arcmin.csv")
    target = rows(recon/"J0754_apass_sources.csv")[0]
    def select(catalog):
        candidates = []
        for row in catalog:
            if not eligible_catalog(row, target):
                continue
            ra, dec = number(row["ra_deg"]), number(row["dec_deg"])
            nearest = min(separation_arcsec(ra, dec, number(other["ra_deg"]), number(other["dec_deg"])) for other in catalog if other["ref_number"] != row["ref_number"])
            if nearest < 20:
                continue
            row = dict(row, nearest_catalog_arcsec=nearest,
                       distance_arcsec=separation_arcsec(ra,dec,118.69125,16.694722222))
            candidates.append(row)
        return sorted(candidates, key=lambda r:r["distance_arcsec"])[:12]
    selected = select(catalog)
    if len(selected) < 8:
        payload = dict(refcat="apass", ra_deg=118.69125, dec_deg=16.694722222, radius_arcsec=1800)
        path = out/"comparison_catalog_30arcmin.csv"
        if not path.exists():
            save_csv(path, post_csv("dasch/dr7/querycat", payload))
        catalog = rows(path)
        selected = select(catalog)
    (out/"selected_controls.json").write_text(json.dumps(selected, indent=2), encoding="utf-8")
    for idx, row in enumerate(selected):
        label = f"control_{idx+1:02}"
        path = out/(label+".csv")
        if not path.exists():
            save_csv(path,post_csv("dasch/dr7/lightcurve",dict(refcat="apass", ref_number=int(row["ref_number"]), gsc_bin_index=int(row["gsc_bin_index"]))))
            time.sleep(.25)
        print(json.dumps(dict(label=label, star=row["ref_text"], rows=len(rows(path)), separation_arcmin=round(row["distance_arcsec"]/60,2))), flush=True)


def applause(out, recon):
    out.mkdir(parents=True,exist_ok=True)
    target = [r for r in rows(recon/"J0754_applause_cone.csv") if r["natmag"]]
    plate_ids = sorted(set(int(r["plate_id"]) for r in target))
    where = "plate_id IN ("+",".join(map(str,plate_ids))+")"
    # Source-table TAP predicates hit the service's statement timeout; its
    # documented spatial-search endpoint provides a complete smaller field.
    url="https://www.plate-archive.org/conesearch/api/applause_dr4.source_calib/?VERB=3&RA=118.69125&DEC=16.694722222&SR=0.166666667"
    field=out/"applause_field_cone10.xml"
    if not field.exists():
        response=requests.get(url,timeout=120)
        response.raise_for_status();field.write_bytes(response.content)
    (out/"applause_cone_request.json").write_text(json.dumps(dict(url=url,radius_arcmin=10,reason="Documented cone-search fallback after TAP source-table timeout; filter desired plates offline."),indent=2),encoding="utf-8")
    queries = {"applause_plates.xml":"SELECT * FROM applause_dr4.plate WHERE "+where,
               "applause_scans.xml":"SELECT * FROM applause_dr4.scan WHERE "+where}
    (out/"applause_metadata_queries.json").write_text(json.dumps(queries,indent=2), encoding="utf-8")
    for name, query in queries.items():
        path = out/name
        if not path.exists():
            response = requests.post("https://www.plate-archive.org/tap/sync",data=dict(REQUEST="doQuery",LANG="ADQL",FORMAT="votable",QUERY=query),timeout=120)
            response.raise_for_status()
            path.write_bytes(response.content)
        print(json.dumps(dict(file=name,bytes=path.stat().st_size)),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["controls","applause"])
    parser.add_argument("--out",type=pathlib.Path,required=True)
    parser.add_argument("--recon",type=pathlib.Path,required=True)
    args = parser.parse_args()
    globals()[args.mode](args.out.resolve(),args.recon.resolve())
