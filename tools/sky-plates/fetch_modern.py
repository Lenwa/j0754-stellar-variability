"""Public modern time-series snapshot for J0754 and the original comparison set."""
import argparse
import csv
import io
import json
import pathlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from astropy.io import fits


def run(previous, recon, out):
    out.mkdir(parents=True,exist_ok=True)
    original=json.loads((previous/"isotonic_pixel_results.json").read_text())
    ids=original["results"][-1]["selected_controls"]+["667401913089980544"]
    with (previous/"gaia_field.csv").open() as f:
        gaia={r["source_id"]:r for r in csv.DictReader(f)}
    def fetch(sid):
        r=gaia[sid]
        path=out/f"ztf_{sid}.csv"
        parameters=dict(CIRCLE=f"{r['ra']} {r['dec']} {3/3600}",FORMAT="csv")
        record=dict(source_id=sid,ra=float(r["ra"]),dec=float(r["dec"]),
                    url="https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves",parameters=parameters)
        try:
            if not path.exists():
                response=requests.get(record["url"],params=parameters,timeout=50)
                response.raise_for_status();path.write_bytes(response.content)
            record["bytes"]=path.stat().st_size
            record["status"]="downloaded_or_cached"
        except requests.RequestException as e:
            record["status"]="request_failed";record["error"]=str(e)
        (out/f"ztf_{sid}_request.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
        return record
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for task in as_completed([pool.submit(fetch,sid) for sid in ids]):
            record=task.result();records.append(record)
            print(record["source_id"],record["status"],record.get("bytes"),flush=True)
    target=next(t for t in json.loads((recon/"recon_manifest.json").read_text())["targets"] if t["label"]=="J0754")
    tess=[]
    for product in target["public_timeseries_metadata"]:
        uri=product["dataURL"]
        if not uri.endswith(".fits"):
            continue
        name=uri.rsplit("/",1)[-1];path=out/name
        record=dict(uri=uri,file=name)
        try:
            if not path.exists():
                response=requests.get("https://mast.stsci.edu/api/v0.1/Download/file",params=dict(uri=uri),timeout=50)
                response.raise_for_status()
                if len(response.content)>30_000_000:
                    raise ValueError("Unexpectedly large lightcurve")
                with fits.open(io.BytesIO(response.content)) as h:
                    assert len(h)>1
                path.write_bytes(response.content)
            record["status"]="downloaded_or_cached";record["bytes"]=path.stat().st_size
        except (requests.RequestException,ValueError,OSError) as e:
            record["status"]="request_failed";record["error"]=str(e)
        tess.append(record);print("TESS",name,record["status"],flush=True)
    (out/"modern_downloads.json").write_text(json.dumps(dict(ztf=sorted(records,key=lambda r:r["source_id"]),tess=tess),indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    for option in ["previous","recon","out"]:p.add_argument("--"+option,type=pathlib.Path,required=True)
    a=p.parse_args();run(a.previous,a.recon,a.out)
