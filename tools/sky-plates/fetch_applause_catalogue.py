"""Retrieve/cache the independent catalogue field and rebuild the small union."""
import argparse
import json
import math
import pathlib
import shutil
import requests
from astropy.io.votable import parse_single_table
from astropy.table import unique, vstack

GAIA_QUERY = ("SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,phot_g_mean_mag,"
              "phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,ruwe FROM gaiadr3.gaia_source "
              "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
              "CIRCLE('ICRS',118.69125,16.694722222,0.32)) AND phot_g_mean_mag<18.5")
UNION_COLUMNS = ["source_id","process_id","scan_id","plate_id","ra_icrs","dec_icrs",
                 "sextractor_flags","nn_dist","model_prediction","cat_natmag","natmag",
                 "gaiaedr3_bp_rp","gaiaedr3_id","gaiaedr3_dist","phot_range_flags","annular_bin"]


def cached_get(session, out, name, url, parameters):
    path = out/name
    if path.exists():
        return path
    r = session.get(url, params=parameters, timeout=55)
    r.raise_for_status()
    path.write_bytes(r.content)
    (out/(path.stem+"_request.json")).write_text(json.dumps(
        dict(url=url,parameters=parameters,bytes=len(r.content)),indent=2)+"\n",encoding="utf-8")
    print(f"Retrieved {name}: {len(r.content)} bytes", flush=True)
    return path


def run(old, out, gaia_radius=.32, gaia_service="esa", registration_from=None):
    out.mkdir(parents=True,exist_ok=True)
    s = requests.Session()
    paths = [old/"applause_field_cone10.xml"]
    if not paths[0].exists():
        raise ValueError("Retained first-pilot 10-arcminute field is required")
    # A single 15-arcminute API query timed out in the saved run. Smaller
    # overlapping fields complete the useful part of the image with less load.
    if registration_from is not None:
        shutil.copyfile(registration_from/"catalogue_union.fits",out/"catalogue_union.fits")
        print("Reused retained registration catalogue",flush=True)
    else:
        for i,(dx,dy) in enumerate([(-6,-6),(-6,6),(6,-6),(6,6)]):
            parameters = dict(RA=118.69125+dx/60/math.cos(math.radians(16.694722222)),
                              DEC=16.694722222+dy/60, SR=.125, VERB=3)
            paths.append(cached_get(s,out,f"catalogue_subfield_{i}.xml",
                "https://www.plate-archive.org/conesearch/api/applause_dr4.source_calib/",parameters))
        merged = unique(vstack([parse_single_table(p).to_table() for p in paths]),keys="source_id")
        merged[UNION_COLUMNS].write(out/"catalogue_union.fits",overwrite=True)
        print(f"Prepared catalogue union: {len(merged)} unique extraction rows",flush=True)
    query=GAIA_QUERY.replace("0.32)",f"{gaia_radius:g})")
    receipt=out/"gaia_field_request.json"
    if (out/"gaia_field.csv").exists() and receipt.exists():
        prior=json.loads(receipt.read_text())
        saved_query=prior.get("query",prior.get("parameters",{}).get("QUERY"))
        if saved_query != query:
            raise ValueError("Cached Gaia field query differs; use a separate output directory")
    url=("https://gea.esac.esa.int/tap-server/tap/sync" if gaia_service=="esa" else
         "https://gaia.ari.uni-heidelberg.de/tap/sync")
    cached_get(s,out,"gaia_field.csv",url,dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query))


if __name__ == "__main__":
    p=argparse.ArgumentParser();p.add_argument("--old",type=pathlib.Path,required=True)
    p.add_argument("--out",type=pathlib.Path,required=True)
    p.add_argument("--gaia-radius",type=float,default=.32)
    p.add_argument("--gaia-service",choices=["esa","heidelberg"],default="esa")
    p.add_argument("--registration-from",type=pathlib.Path)
    a=p.parse_args();run(a.old,a.out,a.gaia_radius,a.gaia_service,a.registration_from)
