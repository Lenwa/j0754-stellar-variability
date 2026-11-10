"""Quality-screened, object-separated modern photometry; no mechanism inference."""
import argparse
import csv
import json
import math
import pathlib
from collections import defaultdict
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def good_ztf(r):
    try:
        return (int(r["catflags"])==0 and 0<float(r["magerr"])<=.1
                and math.isfinite(float(r["mag"])) and math.isfinite(float(r["mjd"])))
    except (ValueError,KeyError):
        return False


def summarize_ztf(rows):
    grouped=defaultdict(list)
    for r in rows:
        if good_ztf(r):grouped[(r["oid"],r["filtercode"])].append(r)
    groups=[]
    for (oid,band),rr in sorted(grouped.items()):
        mag=np.array([float(r["mag"]) for r in rr]);time=np.array([float(r["mjd"]) for r in rr])
        center=float(np.median(mag))
        groups.append(dict(oid=oid,band=band,points=len(rr),median_mag=center,
                           robust_scatter_mag=float(1.4826*np.median(np.abs(mag-center))),
                           p05_p95_mag=np.quantile(mag,[.05,.95]).tolist(),
                           min_max_mag=[float(mag.min()),float(mag.max())],
                           median_error=float(np.median([float(r["magerr"]) for r in rr])),
                           fainter_than_median_by_0_3=int((mag>center+.3).sum()),
                           brighter_than_median_by_0_3=int((mag<center-.3).sum()),
                           first_mjd=float(time.min()),last_mjd=float(time.max()),
                           first_date=Time(time.min(),format="mjd").isot[:10],
                           last_date=Time(time.max(),format="mjd").isot[:10]))
    return dict(raw_points=len(rows),good_points=sum(g["points"] for g in groups),groups=groups)


def run(out):
    ztf={};raw={}
    for path in sorted(out.glob("ztf_*.csv")):
        with path.open() as f:rows=list(csv.DictReader(f))
        sid=path.stem.removeprefix("ztf_")
        if rows and "oid" not in rows[0]:
            ztf[sid]=dict(status="unparsed_response",bytes=path.stat().st_size);continue
        ztf[sid]=dict(status="lightcurve_rows" if rows else "no_rows_in_cone",**summarize_ztf(rows))
        raw[sid]=rows
    tess=[];series=[]
    for path in sorted(out.glob("hlsp_*.fits")):
        with fits.open(path) as h:
            d=h[1].data;cols=d.names
            flux_name=next(name for name in ["PDCSAP_FLUX","KSPSAP_FLUX","DET_FLUX"] if name in cols)
            good=(d["QUALITY"]==0)&np.isfinite(d["TIME"])&np.isfinite(d[flux_name])&(d[flux_name]>0)
            time=np.array(d["TIME"][good],float)+2457000
            flux=np.array(d[flux_name][good],float)
            if not len(flux):continue
            norm=flux/np.median(flux);mag=-2.5*np.log10(norm)
            record=dict(file=path.name,sector=int(h[0].header["SECTOR"]),column=flux_name,
                        raw_points=len(d),good_points=len(flux),
                        first_date=Time(time.min(),format="jd",scale="tdb").isot[:10],
                        last_date=Time(time.max(),format="jd",scale="tdb").isot[:10],
                        robust_scatter_mag=float(1.4826*np.median(np.abs(mag-np.median(mag)))),
                        p05_p95_relative_mag=np.quantile(mag,[.05,.95]).tolist(),
                        crowdsap=h[1].header.get("CROWDSAP"))
            tess.append(record);series.append((record,time,norm))
    fig,axes=plt.subplots(4,1,figsize=(11,10),layout="constrained")
    ids=["667401913089980544","667024467067954560","667023509291889792"]
    labels=["J0754","Comparison absent on 1937 scan (â€¦7954560)","Comparison with negative 1937 aperture (â€¦1889792)"]
    for ax,sid,label in zip(axes[:3],ids,labels):
        rows=raw.get(sid,[])
        for group in ztf.get(sid,{}).get("groups",[]):
            rr=[r for r in rows if good_ztf(r) and r["oid"]==group["oid"] and r["filtercode"]==group["band"]]
            t=Time([float(r["mjd"]) for r in rr],format="mjd").decimalyear
            m=np.array([float(r["mag"]) for r in rr])-group["median_mag"]
            ax.scatter(t,m,s=5,alpha=.6,label=f"{group['band']}, object â€¦{group['oid'][-5:]}")
        ax.axhline(0,color="gray",lw=.7);ax.invert_yaxis();ax.set_title(label,fontsize=11)
        ax.set_ylabel("Relative mag");ax.legend(fontsize=7,loc="best");ax.grid(alpha=.15)
    ax=axes[3]
    # Separate sectors/pipelines by their own median; this cannot measure decade-scale fading.
    for record,time,flux in series:
        if record["column"]=="PDCSAP_FLUX" and record["sector"] in [71,72]:
            xx, yy = time-2460000, (flux-1)*100
            gaps = np.flatnonzero(np.diff(xx) > .1)+1
            xx, yy = np.insert(xx,gaps,np.nan), np.insert(yy,gaps,np.nan)
            ax.plot(xx,yy,lw=.5,label=f"TESS sector {record['sector']}")
    ax.set_xlabel("BJD âˆ’ 2460000");ax.set_ylabel("Relative flux (%)");ax.legend(fontsize=8)
    ax.set_title("J0754 TESS: sector-normalised measurements, 2023",fontsize=11);ax.grid(alpha=.15)
    fig.suptitle("Modern observations: quality-screened, with each object/filter zero point kept separate",fontsize=13)
    fig.savefig(out/"modern_lightcurves.png",dpi=160);plt.close(fig)
    result=dict(ztf=ztf,tess=tess,limits=["No rows in a3arcsec cone is not proof that survey images are absent.",
                "catflags0 and reported errors do not by themselves rule out saturation, blending or all instrumental effects.",
                "Scatter is computed within each object/filter; no cross-detector absolute trend or period claim.",
                "TESS reductions and sector normalisation can suppress long-term variations. Modern stability cannot prove1937stability."])
    (out/"modern_summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for sid in ids:print(sid,json.dumps(ztf.get(sid)))
    print("ZTF sources with rows",sum(v.get("raw_points",0)>0 for v in ztf.values()),"of",len(ztf),"TESS files",len(tess))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True)
    run(p.parse_args().out)
