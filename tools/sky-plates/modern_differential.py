"""Exploratory same-exposure check of modern target variability, preserving tails."""
import argparse
import csv
import json
import pathlib
from collections import defaultdict
import numpy as np
from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from summarize_modern import good_ztf


def frame(r):
    return tuple(r[k] for k in ["expid","filtercode","field","ccdid","qid"])


def run(out):
    sources={}
    for path in out.glob("ztf_*.csv"):
        with path.open() as f:rows=list(csv.DictReader(f))
        sources[path.stem.removeprefix("ztf_")]=[r for r in rows if good_ztf(r)]
    target=sources.pop("667401913089980544")
    centers={}
    for sid,rows in sources.items():
        groups=defaultdict(list)
        for r in rows:groups[(r["oid"],r["filtercode"])].append(float(r["mag"]))
        centers[sid]={k:float(np.median(v)) for k,v in groups.items()}
    controls=defaultdict(dict)
    for sid,rows in sources.items():
        for r in rows:
            key=frame(r);res=float(r["mag"])-centers[sid][(r["oid"],r["filtercode"])]
            # Repeated objects from the same star/frame do not become extra controls.
            value=(float(r["magerr"]),res)
            if sid not in controls[key] or value[0]<controls[key][sid][0]:controls[key][sid]=value
    principal=max({r["oid"] for r in target},key=lambda oid:sum(r["oid"]==oid for r in target))
    center=float(np.median([float(r["mag"]) for r in target if r["oid"]==principal]))
    result=[]
    for r in target:
        if r["oid"]!=principal:continue
        cc=controls.get(frame(r),{})
        if len(cc)<3:continue
        offset=float(np.median([v[1] for v in cc.values()]))
        result.append(dict(mjd=float(r["mjd"]),raw_mag=float(r["mag"]),
                           relative_mag=float(r["mag"])-center,
                           corrected_relative_mag=float(r["mag"])-center-offset,
                           controls=len(cc),control_offset=offset,expid=r["expid"]))
    result.sort(key=lambda r:r["mjd"])
    years=defaultdict(list)
    for r in result:years[int(Time(r["mjd"],format="mjd").decimalyear)].append(r["corrected_relative_mag"])
    annual=[dict(year=y,points=len(v),median_relative_mag=float(np.median(v)),
                 p16_p84=np.quantile(v,[.16,.84]).tolist()) for y,v in sorted(years.items()) if len(v)>=10]
    values=np.array([r["corrected_relative_mag"] for r in result])
    summary=dict(status="Exploratory diagnostic designed after modern scatter/plots were seen; not a pre-specified physical-model test.",
                 principal_object=principal,minimum_controls=3,points=len(result),
                 median_controls=float(np.median([r["controls"] for r in result])),
                 robust_scatter_mag=float(1.4826*np.median(np.abs(values-np.median(values)))),
                 annual=annual,measurements=result,
                 limits=["All available original control stars retained; no control rejected on historical or modern target outcome.",
                         "Image/filter/field/CCD/quadrant key must agree. Same-source duplicate controls collapse to one.",
                         "Detector, saturation, colour and long-term reference-star systematics remain; no physical cause or historical connection inferred.",
                         "Annual quantiles describe observations, not uncertainty of an independent annual mean. Large tails remain and need image vetting."])
    (out/"modern_differential.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    fig,axes=plt.subplots(2,1,figsize=(10,7),layout="constrained")
    tt=Time([r["mjd"] for r in result],format="mjd").decimalyear
    axes[0].scatter(tt,values,s=7,alpha=.6);axes[0].set_title("All accepted measurements, including large tails")
    yy=np.array([r["year"]+.5 for r in annual]);mm=np.array([r["median_relative_mag"] for r in annual])
    error=np.array([[r["median_relative_mag"]-r["p16_p84"][0],r["p16_p84"][1]-r["median_relative_mag"]] for r in annual]).T
    axes[1].errorbar(yy,mm,yerr=error,fmt="o",capsize=3)
    axes[1].set_title("Year medians with observed 16thâ€“84th percentile ranges (10+ points/year)")
    for ax in axes:ax.invert_yaxis();ax.axhline(0,color="gray",lw=.7);ax.grid(alpha=.2);ax.set_ylabel("Relative g magnitude")
    axes[1].set_xlabel("Year")
    fig.suptitle("J0754: same-exposure comparison-star correction in the largest ZTF data set",fontsize=12)
    fig.savefig(out/"modern_differential.png",dpi=160);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k not in ["measurements","limits"]},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True)
    run(p.parse_args().out)
