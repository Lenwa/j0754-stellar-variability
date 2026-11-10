"""Public IRSA snapshots for J0754 and its original fixed comparison set."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import datetime
import json
from pathlib import Path

import requests
from astropy.table import Table

TARGET = "667401913089980544"


def gator(out, name, table, ra, dec, radius, columns="*"):
    """Preserve native IPAC output, including column types and leading-zero flags."""
    native=out/(name+".tbl")
    if not native.exists():
        params=dict(catalog=table,spatial="Cone",objstr=f"{ra} {dec}",
                    radius=radius,radunits="arcsec",outfmt=1)
        if columns!="*":
            params["selcols"]=columns
        response=requests.get("https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query",
                              params=params,timeout=55)
        (out/(name+"_gator_request.json")).write_text(json.dumps(dict(
            params=params,url=response.url,status=response.status_code,
            retrieved_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            bytes=len(response.content)),indent=2)+"\n")
        response.raise_for_status()
        # Parse before accepting a service error page as a scientific table.
        Table.read(response.text.splitlines(),format="ascii.ipac")
        native.write_bytes(response.content)
    data=Table.read(native,format="ascii.ipac")
    if len(data)>=100000:
        raise ValueError("Possible truncated table: "+name)
    data.write(out/(name+".csv"),format="ascii.csv",overwrite=True)
    return dict(file=name+".csv",native=native.name,rows=len(data),bytes=native.stat().st_size)


def run(previous,out):
    out.mkdir(exist_ok=True,parents=True)
    ids=json.loads((previous/"isotonic_pixel_results.json").read_text())["results"][-1]["selected_controls"]+[TARGET]
    with (previous/"gaia_field.csv").open() as f:
        gaia={r["source_id"]:r for r in csv.DictReader(f)}
    columns="ra,dec,source_id,cntr,mjd,scan_id,frame_num,w1mpro,w1sigmpro,w1snr,w1rchi2,w2mpro,w2sigmpro,w2snr,w2rchi2,w1sat,w2sat,nb,na,cc_flags,ph_qual,sso_flg,qual_frame,qi_fact,saa_sep,moon_masked,w1frtr,w2frtr,allwise_cntr,r_allwise,n_allwise,w1mpro_allwise,w2mpro_allwise,w3mpro_allwise,w4mpro_allwise,w3sigmpro_allwise,w4sigmpro_allwise"
    def one(sid):
        r=gaia[sid];ra=float(r["ra"]);dec=float(r["dec"]);radius=10 if sid==TARGET else 3
        result=gator(out,"neowise_"+sid,"neowiser_p1bs_psd",ra,dec,radius,columns)
        result.update(gaia_id=sid,ra=ra,dec=dec,radius_arcsec=radius)
        print(sid,result["rows"],"rows",flush=True)
        return result
    records=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(one,sid) for sid in ids]):
            records.append(future.result())
    r=gaia[TARGET];ra=float(r["ra"]);dec=float(r["dec"])
    extra=[]
    for table,name,radius in [("allwise_p3as_psd","allwise_target",10),
                               ("allwise_p3as_mep","allwise_epochs_target",3),
                               ("fp_psc","twomass_target",3)]:
        result=gator(out,name,table,ra,dec,radius);extra.append(result)
        print(name,result["rows"],"rows",flush=True)
    (out/"infrared_downloads.json").write_text(json.dumps(dict(neowise=sorted(records,key=lambda r:r["gaia_id"]),other=extra),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--previous",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.previous,a.out)
