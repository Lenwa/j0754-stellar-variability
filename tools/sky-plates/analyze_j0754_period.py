"""Test the catalogue clock against later data; longer-period scans are exploratory."""
import argparse
import csv
import json
from pathlib import Path

from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from asassn_replication import parse_row
from period_models import PERIOD, EPOCH, bin_series, fit_model, design, scores, forecast, blocked_cv, replay_pairs

CUT = float(Time("2023-01-01").jd)
START = float(Time("2019-01-01").jd)


def arrays(bins):
    return (np.array([r["jd"] for r in bins]),np.array([r["mag"] for r in bins]),
            np.array([r["camera"] for r in bins]))


def basic(bins):
    t,y,c=arrays(bins)
    null=fit_model(t,y,c,"null")["prediction"]
    result=dict(bins=len(t),cameras=sorted(set(c)),baseline_days=float(np.ptp(t)),models={})
    for label,kind,harmonics in [("one_harmonic","periodic",1),("two_harmonics","periodic",2),("slow_reference","slow",1)]:
        model=fit_model(t,y,c,kind,harmonics=harmonics)
        result["models"][label]=dict(full=scores(y,model["prediction"],null),parameters=model["parameters"],
            blocked_cv=blocked_cv(t,y,c,kind,harmonics=harmonics),
            forward_forecast=(forecast(t,y,c,kind,CUT,harmonics=harmonics) if kind=="periodic" else
                dict(status="NOT_IDENTIFIABLE_FOR_EXTRAPOLATION",reason="Future knot values are not identified by pre2023training; slow curve is a descriptive/interpolation reference only.")))
    pairs=replay_pairs(t,y,c)
    values=np.array([r["delta"] for r in pairs])
    result["two_cycle_replay"]=dict(pairs=pairs,count=len(pairs),
        rms=float(np.sqrt(np.mean(values**2))) if len(values) else None,
        median_absolute=float(np.median(np.abs(values))) if len(values) else None,
        by_camera={str(cam):dict(count=len(v),rms=float(np.sqrt(np.mean(v**2))))
                   for cam in sorted(set(c)) if len(v:=np.array([r["delta"] for r in pairs if r["camera"]==cam]))})
    result["per_camera_two_harmonic_r2"]={}
    for cam in sorted(set(c)):
        use=c==cam
        local=fit_model(t[use],y[use],c[use],"periodic",harmonics=2)
        baseline=fit_model(t[use],y[use],c[use],"null")
        result["per_camera_two_harmonic_r2"][str(cam)]=scores(y[use],local["prediction"],baseline["prediction"])["r2"]
    return result


def extra(bins,seed):
    t,y,c=arrays(bins);null=fit_model(t,y,c,"null")["prediction"]
    frequency=np.linspace(1/6000,1/200,4096);periods=1/frequency
    power=np.array([scores(y,fit_model(t,y,c,"periodic",period=p)["prediction"],null)["r2"] for p in periods])
    window=np.abs(np.mean(np.exp(2j*np.pi*frequency[:,None]*(t-EPOCH)),axis=1))**2
    best=int(np.argmax(power))
    train=t<CUT
    nearby=np.linspace(.9*PERIOD,1.1*PERIOD,401)
    training_error=[np.mean((y[train]-fit_model(t,y,c,"periodic",period=p,harmonics=2,train=train)["prediction"][train])**2) for p in nearby]
    local=float(nearby[int(np.argmin(training_error))])
    # Stable-signal power experiment at the actual retained time/camera bins.
    rng=np.random.default_rng(seed);blocks=np.floor((t-EPOCH)/90).astype(int);unique,inv=np.unique(blocks,return_inverse=True)
    offset=np.array([12+.05*sorted(set(c)).index(cam) for cam in c])
    successful=0;simulated_replays=[];gains=[];variances=[]
    for _ in range(200):
        phase=rng.uniform(0,2*np.pi)
        signal=offset+.15*np.sin(2*np.pi*(t-EPOCH)/PERIOD+phase)
        synthetic=signal+rng.normal(0,.03,len(unique))[inv]+rng.normal(0,.04,len(t))
        fit=fit_model(t,synthetic,c,"periodic")
        base=fit_model(t,synthetic,c,"null")
        r2=scores(synthetic,fit["prediction"],base["prediction"])["r2"]
        gain=forecast(t,synthetic,c,"periodic",CUT)["rmse_gain"]
        successful += r2>=.5 and gain>=.2;gains.append(gain);variances.append(r2)
        pairs=replay_pairs(t,synthetic,c)
        simulated_replays.append(float(np.sqrt(np.mean([r["delta"]**2 for r in pairs]))))
    return dict(exploratory_scan=dict(period_days=periods.tolist(),camera_centered_r2=power.tolist(),window_power=window.tolist(),
                best_period=float(periods[best]),best_r2=float(power[best]),baseline_cycles=float(np.ptp(t)/periods[best]),
                boundary=best in (0,len(periods)-1),warning="Exploratory profile only; no false-alarm probability or confirmed replacement period."),
                nearby_training_selected=dict(period=local,forecast=forecast(t,y,c,"periodic",CUT,period=local,harmonics=2)),
                cadence_simulation=dict(trials=200,successful=int(successful),fraction=float(successful/200),
                    median_r2=float(np.median(variances)),median_forecast_rmse_gain=float(np.median(gains)),
                    replay_rms_p95=float(np.quantile(simulated_replays,.95)),
                    assumptions="0.30mag peak-to-peak stable sinusoid, 0.03mag shared90-day noise plus0.04mag independent bin noise; conditional power, not observed false-alarm rate."))


def decide(results):
    supported=[];failed=[]
    for r in results.values():
        models=[r["models"][name] for name in ("one_harmonic","two_harmonics")]
        supported.append(any(m["full"]["r2"]>=.5 and m["blocked_cv"]["rmse_gain"]>=.2
                             and m["forward_forecast"]["rmse_gain"]>=.2 for m in models))
        weak=all(m["full"]["r2"]<.5 and m["blocked_cv"]["rmse_gain"]<.2
                 and m["forward_forecast"]["rmse_gain"]<.2 for m in models)
        weak &= r["nearby_training_selected"]["forecast"]["rmse_gain"]<.2
        failed.append(weak and r["two_cycle_replay"]["rms"]>max(.1,r["cadence_simulation"]["replay_rms_p95"]))
    return "SUPPORTED_DOMINANT_912_CYCLE" if all(supported) else "STABLE_DOMINANT_912_MODEL_FAILS" if all(failed) else "MIXED_OR_INCONCLUSIVE"


def run(independent,ztf,out):
    out.mkdir(parents=True,exist_ok=True)
    with (independent/"asassn_v1_lightcurve.csv").open(newline="") as f:
        all_rows=[parse_row(r) for r in csv.DictReader(f)]
    g=[r for r in all_rows if r["accepted"] and r["band"]=="g" and r["jd"]>=START]
    z=json.loads((ztf/"modern_differential.json").read_text())
    zrows=[]
    for r in z["measurements"]:
        jd=r["mjd"]+2400000.5
        if jd>=START:
            zrows.append(dict(jd=jd,mag=r["corrected_relative_mag"],day=Time(jd,format="jd").iso[:10],
                              camera="ZTF_"+z["principal_object"],band="g"))
    datasets={"ASAS_SN_g":g,"ZTF_g":zrows}
    results={};binned={}
    for i,(name,rows) in enumerate(datasets.items()):
        bins=bin_series(rows);binned[name]=bins
        result=basic(bins);result.update(extra(bins,20261108+i))
        result["native_accepted_points"]=len(rows)
        result["bin_width_sensitivities"]={str(width):basic(bin_series(rows,width)) for width in (15,60)}
        if name=="ASAS_SN_g":
            strict=[r for r in rows if r["mag_err"]<=.1 and 0<r["fwhm"]<=2.5]
            result["strict_quality_sensitivity"]=basic(bin_series(strict))
        results[name]=result
        print(name,"bins",len(bins),"two-harmonic R2",result["models"]["two_harmonics"]["full"]["r2"],flush=True)
    old=json.loads((independent/"asassn_variable_target.json").read_text())["results"]
    legacy=[dict(jd=r["hjd"],mag=r["mag"],camera=r["camera"],band="V",day=Time(r["hjd"],format="jd").iso[:10])
            for r in old if 0<r["mag_err"]<99 and r["flux"]>0 and r["flux_err"]>0]
    t,y,c=arrays(bin_series(legacy));base=fit_model(t,y,c,"null")
    legacy_result=dict(native_points=len(legacy),bins=len(t),note="Original V fitting data, not independent validation.",
        models={str(h):scores(y,fit_model(t,y,c,"periodic",harmonics=h)["prediction"],base["prediction"]) for h in (1,2)})
    verdict=decide(results)
    summary=dict(period_days=PERIOD,epoch_hjd=EPOCH,verdict=verdict,datasets=results,legacy_V=legacy_result,
        limits=["Broad modern trend was already known; protocol predates model fitting, not all data inspection.",
                "Tests address a stable dominant cycle, not absence of a weak/quasiperiodic component.",
                "Binning/camera offsets and time-block validation reduce common pitfalls but do not remove all calibration/correlated errors.",
                "HJD versus geocentric ZTF JD difference is negligible for these month-scale checks; no precision timing/orbit inference.",
                "Longer profile peaks over few cycles do not measure a replacement period or establish rotation, dust or a companion."])
    (out/"period_results.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    for name,bins in binned.items():
        with (out/(name+"_bins.csv")).open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(bins[0]));writer.writeheader();writer.writerows(bins)
    sensitivity={}
    for name,bins in binned.items():
        t,y,c=arrays(bins);sensitivity[name]={}
        for tolerance in (5,10,15):
            pairs=replay_pairs(t,y,c,tolerance=tolerance);delta=np.array([r["delta"] for r in pairs])
            sensitivity[name][str(tolerance)]=dict(pairs=len(pairs),
                rms=float(np.sqrt(np.mean(delta**2))) if len(delta) else None,
                median_absolute=float(np.median(np.abs(delta))) if len(delta) else None)
    (out/"replay_tolerance_sensitivity.json").write_text(json.dumps(dict(
        status="Post-result check of finite phase-matching tolerance; primary15-day tolerance unchanged.",results=sensitivity),indent=2)+"\n")
    plot(summary,binned,out)
    print("Verdict",verdict,flush=True)


def plot(summary,binned,out):
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout="constrained")
    bins=[r for r in binned["ASAS_SN_g"] if r["camera"]=="bi"]
    t,y,c=arrays(bins);grid=np.linspace(t.min(),t.max(),1200);gc=np.array(["bi"]*len(grid))
    periodic=fit_model(t,y,c,"periodic",harmonics=2)
    slow=fit_model(t,y,c,"slow")
    axes[0,0].scatter(Time(t,format="jd").decimalyear,y,s=22,color="tab:blue",label="ASAS-SN bi, 30-day bins")
    axes[0,0].plot(Time(grid,format="jd").decimalyear,design(grid,gc,"periodic",harmonics=2)@periodic["beta"],color="tab:red",label="Best fixed 912-day fit (2 harmonics)")
    axes[0,0].plot(Time(grid,format="jd").decimalyear,design(grid,gc,"slow",knots=np.linspace(t.min(),t.max(),5))@slow["beta"],"--",color="grey",label="Slow reference, same parameter count")
    axes[0,0].set(title="The observed broad change versus a repeating clock",xlabel="Year (2026 is partial)",ylabel="g magnitude")
    axes[0,0].invert_yaxis();axes[0,0].legend(fontsize=8)
    cycles=np.floor((t-EPOCH)/PERIOD).astype(int);phase=((t-EPOCH)/PERIOD)%1
    pts=axes[0,1].scatter(phase,y,c=cycles,cmap="viridis",s=28)
    fig.colorbar(pts,ax=axes[0,1],label="Cycle index from catalogue epoch",ticks=sorted(set(cycles)))
    axes[0,1].set(title="Same phases have different brightness in different cycles",xlabel="Phase at 912.3159 days",ylabel="g magnitude (one camera, bi)")
    axes[0,1].invert_yaxis()
    for name,r in summary["datasets"].items():
        scan=r["exploratory_scan"]
        axes[1,0].plot(scan["period_days"],scan["camera_centered_r2"],label=name.replace("ASAS_SN","ASAS-SN").replace("_"," "))
    scan=summary["datasets"]["ASAS_SN_g"]["exploratory_scan"]
    axes[1,0].plot(scan["period_days"],scan["window_power"],color="grey",alpha=.35,label="ASAS-SN sampling window")
    axes[1,0].axvline(PERIOD,color="red",ls="--",label="Catalogue: 912 days")
    axes[1,0].set(xscale="log",xlabel="Trial period (days)",ylabel="Variance explained / window power",title="Exploratory single-harmonic scan",ylim=(-.02,1))
    axes[1,0].minorticks_off()
    axes[1,0].set_xticks([200,365.25,912,2000,4000,6000],labels=["200","365","912","2000","4000","6000"])
    axes[1,0].legend(fontsize=8)
    for name,r in summary["datasets"].items():
        pairs=r["two_cycle_replay"]["pairs"]
        axes[1,1].scatter(Time([p["first_jd"] for p in pairs],format="jd").decimalyear,[p["delta"] for p in pairs],s=20,alpha=.65,label=name.replace("ASAS_SN","ASAS-SN").replace("_"," "))
    axes[1,1].axhline(0,color="black",lw=.8)
    axes[1,1].axhspan(-.1,.1,color="grey",alpha=.1,label="Â±0.10mag reference")
    axes[1,1].set(title="Brightness after two cycles (about five years)",xlabel="Earlier observation year",ylabel="Later minus earlier magnitude (same camera)")
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle("J0754: testing the catalogue's 912-day cycle",fontsize=14)
    fig.savefig(out/"period_test.png",dpi=160);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("independent","ztf","out"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();run(a.independent,a.ztf,a.out)
