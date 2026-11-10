"""Retrieve and render same-series before/event/recovery scans for defect checks."""
import argparse
import base64
import csv
import json
import pathlib
import requests
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run(out):
    with (out/"dasch_differential.csv").open() as stream: rows=list(csv.DictReader(stream))
    rows=[r for r in rows if r["series"]=="rh"]
    selections=[]
    for name,lo,hi,reference in [("Before",1910,1930,1928),("Faded",1930,1936,1933),("Recovery",1941,1945,1942)]:
        row=min((r for r in rows if lo<=float(r["year"])<hi),key=lambda r:abs(float(r["year"])-reference))
        _,plate,_,sol=row["key"].split(":")
        payload=dict(plate_id="rh"+plate.zfill(5),solution_number=int(sol),center_ra_deg=118.69125,center_dec_deg=16.694722222)
        path=out/("scan_"+name.lower()+".fits.gz")
        if not path.exists():
            response=requests.post("https://api.starglass.cfa.harvard.edu/public/dasch/dr7/cutout",json=payload,timeout=60)
            response.raise_for_status();path.write_bytes(base64.b64decode(response.json()))
        selections.append(dict(name=name,path=str(path),payload=payload,year=float(row["year"]),B=float(row["mag"])))
        print(json.dumps(selections[-1]),flush=True)
    selected=json.loads((out/"selected_controls.json").read_text())[0]
    fig,axes=plt.subplots(1,3,figsize=(13,4.8),layout="constrained")
    for ax,item in zip(axes,selections):
        with fits.open(item["path"]) as hdus:
            data=hdus[0].data.astype(float);wcs=WCS(hdus[0].header)
        x,y=wcs.world_to_pixel_values(118.69125,16.694722222)
        cx,cy=wcs.world_to_pixel_values(float(selected["ra_deg"]),float(selected["dec_deg"]))
        norm=ImageNormalize(data,interval=PercentileInterval(99),stretch=AsinhStretch())
        ax.imshow(data,origin="lower",cmap="gray",norm=norm)
        ax.add_patch(plt.Circle((x,y),12,fill=False,color="#e44b41",lw=1.5))
        ax.text(x+16,y+15,"J0754",color="#e44b41",fontsize=10,fontweight="bold")
        ax.add_patch(plt.Circle((cx,cy),12,fill=False,color="#16bfd1",lw=1.5))
        ax.text(cx+16,cy+15,"Control 1",color="#16bfd1",fontsize=10,fontweight="bold")
        ax.set_xlim(x-200,x+200);ax.set_ylim(y-200,y+200);ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(f'{item["name"]}: {item["year"]:.1f}\n{item["payload"]["plate_id"]}, catalog B={item["B"]:.2f}')
    fig.suptitle("Representative RH-series scans â€” displays scaled separately for defect checks",fontsize=13)
    fig.savefig(out/"j0754_scan_check.png",dpi=180);plt.close(fig)
    (out/"scan_selection.json").write_text(json.dumps(selections,indent=2),encoding="utf-8")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=pathlib.Path,required=True)
    run(parser.parse_args().out.resolve())
