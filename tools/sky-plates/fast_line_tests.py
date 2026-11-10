"""Retrospective FAST line robustness tests; counts diagnostics, no atmospheric/RV fit.

Run with --input DIR --out DIR. Inputs are untouched extracted IRAF FITS.
Native-pixel sigma assumes independent noise; convolution/resampling covariance
and fitted-continuum covariance are propagated exactly to first order.
No calibrated intrinsic-variability or temperature inference is supplied.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

FWHM_SIGMA = 2.3548200450309493


def bin_weights(x, interval):
    x = np.asarray(x, float)
    edges = np.r_[x[0]-(x[1]-x[0])/2, (x[:-1]+x[1:])/2,
                  x[-1]+(x[-1]-x[-2])/2]
    return np.maximum(0, np.minimum(edges[1:], interval[1])
                      - np.maximum(edges[:-1], interval[0]))


def equivalent_width(x, y, covariance, line, sides, degree=1, continuum_bias=0):
    """EW and delta-method variance, including continuum/line cross-covariance."""
    x, y, C = np.asarray(x,float), np.asarray(y,float), np.asarray(covariance,float)
    w = bin_weights(x, line)
    lm = w > 0
    cm = np.zeros(len(x),bool)
    for a,b in sides:
        cm |= (x >= a) & (x <= b)
    if np.any(lm & cm):
        raise ValueError("Line and continuum pixel windows overlap")
    if cm.sum() < degree+2 or not lm.any():
        raise ValueError("Insufficient line/continuum samples")
    center = np.mean(line)
    X = np.polynomial.polynomial.polyvander((x-center)/50, degree)
    weights = 1/np.maximum(np.diag(C)[cm], np.finfo(float).tiny)
    Xc = X[cm]
    Bc = np.linalg.solve(Xc.T@(weights[:,None]*Xc), Xc.T*weights)
    B = np.zeros((degree+1,len(x)))
    B[:,cm] = Bc
    beta = B@y
    continuum = X@beta*(1+continuum_bias)
    if np.any(continuum[lm] <= 0):
        raise ValueError("Non-positive continuum")
    ew = np.sum(w[lm]*(1-y[lm]/continuum[lm]))
    gl = np.zeros(len(x)); gl[lm] = -w[lm]/continuum[lm]
    # B fits the un-biased coefficients, hence the explicit multiplicative factor.
    gbeta = np.sum((w[lm]*y[lm]/continuum[lm]**2)[:,None]
                   * X[lm]*(1+continuum_bias),axis=0)
    gc = gbeta@B
    gradient = gl+gc
    vl = float(gl@C@gl)
    vc = float(gc@C@gc)
    cross = float(2*gl@C@gc)
    variance = max(0,float(gradient@C@gradient))
    return {"ew":float(ew),"sigma":float(np.sqrt(variance)),
            "variance_line":vl,"variance_continuum":vc,"variance_cross":cross,
            "gradient":gradient,"continuum":continuum,
            "beta_covariance":B@C@B.T,
            "window_A":float(w.sum())}


def resolution_operator(x, out, native_fwhm, target_fwhm):
    """Linear Gaussian degradation, or interpolation for zero added width."""
    x, out = np.asarray(x,float), np.asarray(out,float)
    if target_fwhm < native_fwhm:
        raise ValueError("Cannot deconvolve to a narrower target")
    sigma = np.sqrt(max(0,target_fwhm**2-native_fwhm**2))/FWHM_SIGMA
    if sigma < 1e-8:
        A = np.zeros((len(out),len(x)))
        for j,t in enumerate(out):
            k = int(np.clip(np.searchsorted(x,t),1,len(x)-1))
            f = (t-x[k-1])/(x[k]-x[k-1])
            A[j,k-1:k+1] = [1-f,f]
        return A
    distance = out[:,None]-x[None,:]
    A = np.exp(-.5*(distance/sigma)**2)*np.gradient(x)[None,:]
    A[np.abs(distance)>6*sigma]=0
    norm=A.sum(axis=1)
    if np.any(norm==0):
        raise ValueError("Output outside kernel support")
    return A/norm[:,None]


def apply_operator(A, y, variance):
    return A@y, (A*np.asarray(variance)[None,:])@A.T


def night_summary(records):
    result={}
    for night in sorted({r["night"] for r in records}):
        rows=[r for r in records if r["night"]==night]
        n=len(rows)
        result[night]={"ew":float(np.mean([r["ew"] for r in rows])),
                      "sigma_formal":float(np.sqrt(sum(r["sigma"]**2 for r in rows))/n),
                      "n_exposures":n,
                      "repeat_sample_sd":float(np.std([r["ew"] for r in rows],ddof=1)) if n>1 else None}
    return result


def independent_night_bounds(records, baseline=("2008-10-22","2008-10-24"),
                             event="2009-11-12"):
    """Conservative extrema over per-exposure scenarios, averaged per night.

    This is a sensitivity envelope, not a statistical confidence interval.
    Repeated exposures never become separate baseline/event nights.
    """
    bounds={}
    for night in (*baseline,event):
        rows=[r for r in records if r["night"]==night]
        files={r["file"] for r in rows}
        if not files:
            raise ValueError("Missing required night")
        bounds[night]=[
            float(np.mean([min(r["ew"] for r in rows if r["file"]==f) for f in files])),
            float(np.mean([max(r["ew"] for r in rows if r["file"]==f) for f in files]))]
    return {"delta_min_A":float(np.mean([bounds[n][0] for n in baseline])-bounds[event][1]),
            "delta_max_A":float(np.mean([bounds[n][1] for n in baseline])-bounds[event][0]),
            "night_bounds":bounds}


def load_spectra(folder):
    from astropy.io import fits
    from astropy.time import Time
    records=[]
    for p in sorted(Path(folder).glob("*.fits")):
        with fits.open(p) as hdus:
            h=hdus[0].header.copy();data=np.array(hdus[0].data,float)
        if data.shape!=(4,1,2681) or not str(h["BANDID4"]).startswith("sigma -"):
            raise ValueError("Unverified IRAF sigma/data bands")
        if h["DC-FLAG"]!=0 or h["WAT2_001"]!="wtype=linear" or "angstroms" not in h["WAT1_001"]:
            raise ValueError("Unsupported wavelength solution")
        x=h["CRVAL1"]+(np.arange(data.shape[-1])+1-h["CRPIX1"])*h["CD1_1"]
        date=Time(float(h["MJD"]),format="mjd",scale="utc").isot
        records.append({"file":p.name,"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
                        "date":date,"night":date[:10],"grating":str(h["DISPERSE"]).strip(),
                        "slit_arcsec":float(h["APERTURE"]),"exptime":float(h["EXPTIME"]),
                        "x":x,"y":data[0,0],"sky":data[2,0],"var":data[3,0]**2,
                        "airmass":str(h.get("AIRMASS"))})
    return sorted(records,key=lambda r:r["date"])


def setup_region(record, bounds, target=10, native=None, offset=0):
    r=record
    native = native if native is not None else (7 if r["grating"]=="300" else 3)
    m=(r["x"]>=bounds[0]-40)&(r["x"]<=bounds[1]+40)
    x=r["x"][m]+offset
    out=np.arange(bounds[0],bounds[1]+.1,.5)
    A=resolution_operator(x,out,native,target)
    y,C=apply_operator(A,r["y"][m],r["var"][m])
    return out,y,C,A@r["sky"][m],A,r["y"][m],r["var"][m]


def scalar_result(r):
    return {k:v for k,v in r.items() if k not in ("gradient","continuum","beta_covariance")}


def contrast(rows):
    nights=night_summary(rows)
    baseline=[nights[x] for x in ["2008-10-22","2008-10-24"]]
    event=nights["2009-11-12"]
    d=float(np.mean([r["ew"] for r in baseline])-event["ew"])
    sigma=float(np.sqrt(sum(r["sigma_formal"]**2 for r in baseline)/4+event["sigma_formal"]**2))
    return {"delta_2008minus2009_A":d,"sigma_formal_A":sigma,"formal_ratio_not_discovery_sigma":d/sigma,
            "nights":nights}


def run(input_dir, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    if not (out/"protocol.txt").exists():
        raise ValueError("Save protocol before real-data comparisons")
    records=load_spectra(input_dir)
    sidesets=[[(6520,6540),(6580,6600)],[(6510,6540),(6580,6610)],
              [(6525,6545),(6590,6610)]]
    primary=[];sens=[];systematics=[];broad=[];curves=[];mc=[]
    rng=np.random.default_rng(20261108)
    for r in records:
        x,y,C,sky,A,yn,vn=setup_region(r,(6500,6620))
        line=(6547.8,6577.8)
        measured=equivalent_width(x,y,C,line,sidesets[0])
        row={"file":r["file"],"night":r["night"],"date":r["date"],"grating":r["grating"],
             **scalar_result(measured)}
        primary.append(row);curves.append((r,x,y/measured["continuum"],sky/measured["continuum"]))
        np.savez_compressed(out/(r["file"]+".common.npz"),wavelength_A=x,signal=y,
                            covariance=C,sky=sky,continuum=measured["continuum"])
        # Validate the real-data sigma at two useful S/N levels, with native-noise injections.
        if r["file"].startswith(("0123.","0146.")):
            draws=(yn+rng.normal(size=(800,len(yn)))*np.sqrt(vn))@A.T
            ws=np.array([equivalent_width(x,z,C,line,sidesets[0])["ew"] for z in draws])
            mc.append({"file":r["file"],"draws":len(ws),"predicted_sigma_A":measured["sigma"],
                       "empirical_sigma_A":float(ws.std(ddof=1)),
                       "mean_minus_measured_A":float(ws.mean()-measured["ew"])})
        for label,values in [("broadband",[-.03,-.015,.015,.03]),
                             ("continuum_bias",[-.03,-.015,.015,.03]),
                             ("sky_band",[-1,-.25,-.1,-.05,.05,.1,.25,1]),
                             ("line_emission",[.1,.3,.5])]:
            for value in values:
                z=y.copy();bias=0
                if label=="broadband":z+=value*measured["continuum"]
                elif label=="continuum_bias":bias=value
                elif label=="sky_band":z+=value*sky
                else:
                    profile=np.exp(-.5*((x-6562.8)/(10/FWHM_SIGMA))**2)
                    profile/=np.sum(profile)*.5
                    z+=value*np.median(measured["continuum"])*profile
                result=equivalent_width(x,z,C,line,sidesets[0],continuum_bias=bias)
                systematics.append({"file":r["file"],"night":r["night"],"type":label,
                                    "amount":value,**scalar_result(result),
                                    "delta_EW_A":result["ew"]-measured["ew"]})
        native_grid=[5,6,7,8] if r["grating"]=="300" else [2,2.3,3,4]
        # Resolution sensitivity uses one primary continuum/window; continuum sensitivity
        # uses primary 10A/7-or3A degradation, avoiding conflating the two.
        for target in [8,10,12]:
            for native in native_grid:
                xx,yy,CC,_,_,_,_=setup_region(r,(6500,6620),target,native)
                result=equivalent_width(xx,yy,CC,line,sidesets[0])
                sens.append({"file":r["file"],"night":r["night"],"family":"LSF",
                             "target":target,"native":native,**scalar_result(result)})
        for half in [10,15,20]:
            for j,sides in enumerate(sidesets):
                for degree in [0,1,2]:
                    try:result=equivalent_width(x,y,C,(6562.8-half,6562.8+half),sides,degree)
                    except ValueError:continue
                    sens.append({"file":r["file"],"night":r["night"],"family":"continuum",
                                 "half":half,"sideset":j,"degree":degree,**scalar_result(result)})
        for shift in [-2,2]:
            xx,yy,CC,_,_,_,_=setup_region(r,(6500,6620),offset=shift)
            result=equivalent_width(xx,yy,CC,line,sidesets[0])
            sens.append({"file":r["file"],"night":r["night"],"family":"wavelength",
                         "offset":shift,**scalar_result(result)})
        for name,iline,sides,bounds in [
            ("Mg_b",(5160,5190),[(5130,5150),(5200,5220)],(5115,5230)),
            ("Na_D",(5880,5906),[(5850,5870),(5910,5930)],(5840,5940))]:
            xx,yy,CC,_,_,_,_=setup_region(r,bounds)
            result=equivalent_width(xx,yy,CC,iline,sides)
            broad.append({"file":r["file"],"night":r["night"],"index":name,**scalar_result(result)})
    base=contrast(primary)
    contrasts=[]
    keys=sorted({tuple((k,str(v)) for k,v in s.items() if k in
                        ("family","target","native","half","sideset","degree","offset")) for s in sens})
    # Native-resolution grids differ between gratings; primary contrast is all300
    # so each grid entry is well defined for the comparison nights.
    for key in keys:
        rows=[s for s in sens if all(str(s.get(k))==v for k,v in key)]
        if {"2008-10-22","2008-10-24","2009-11-12"}.issubset({s["night"] for s in rows}):
            contrasts.append({"specification":dict(key),**contrast(rows)})
    b_ew=np.mean([base["nights"][n]["ew"] for n in ["2008-10-22","2008-10-24"]])
    e_ew=base["nights"]["2009-11-12"]["ew"]
    summary={"primary":base,"continuum_only_beta_to_raise_2009_to_2008":(b_ew-e_ew)/(30-b_ew),
             "broadband_alpha_to_dilute_2008_to_2009":b_ew/e_ew-1,
             "source_files":[{k:r[k] for k in ("file","sha256","date","night","grating","slit_arcsec","exptime","airmass")} for r in records],
             "real_noise_MC":mc,"limitations":["Formal errors assume independent original sigma pixels; unknown native resampling covariance remains.",
              "Per-exposure LSFs are assumptions and diagnostic sensitivity envelopes are not confidence intervals.",
              "Exported background-band amplitude/normalization is not certified residual contamination.",
              "No calibrated continuum temperatures, resolved LDRs or causal/activity determination."]}
    robustness={}
    for family in ["LSF","continuum","wavelength"]:
        variants=[c for c in contrasts if c["specification"]["family"]==family]
        robustness[family]={"variant_count":len(variants),
            "minimum_contrast_A":min(c["delta_2008minus2009_A"] for c in variants),
            "maximum_contrast_A":max(c["delta_2008minus2009_A"] for c in variants)}
    robustness["independent_LSF_by_target"]={str(t):independent_night_bounds(
        [r for r in sens if r["family"]=="LSF" and r["target"]==t]) for t in [8,10,12]}
    robustness["independent_systematic_scenarios"]={name:independent_night_bounds(
        [r for r in systematics if r["type"]==name]) for name in
        ["broadband","continuum_bias","sky_band","line_emission"]}
    robustness["broad_indices"]={name:contrast([r for r in broad if r["index"]==name])
                                 for name in ["Mg_b","Na_D"]}
    for name,value in [("primary",primary),("sensitivities",sens),("systematic_injections",systematics),
                       ("broad_indices",broad),("contrasts",contrasts),("summary",summary),
                       ("robustness_summary",robustness)]:
        (out/(name+".json")).write_text(json.dumps(value,indent=2),encoding="utf-8")
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True)
    for r,x,norm,sky in curves:
        if r["grating"]=="300":
            axs[0].plot(x,norm,label=r["date"][:19],lw=.8)
            axs[1].plot(x,sky,label=r["date"][:10],lw=.8)
    axs[0].set_ylabel("Source / fitted continuum");axs[0].legend(fontsize=7,ncol=2)
    axs[1].set_ylabel("Exported sky / source continuum");axs[1].set_xlabel("Wavelength (Angstrom)")
    for ax in axs:ax.set_xlim(6510,6610);ax.grid(alpha=.2)
    fig.suptitle("FAST300: assumed common10A; exported sky is not a measured subtraction residual")
    fig.tight_layout();fig.savefig(out/"source_sky.png",dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5))
    nights=base["nights"];labels=list(nights)
    ax.errorbar(range(len(labels)),[nights[k]["ew"] for k in labels],
                yerr=[nights[k]["sigma_formal"] for k in labels],fmt="o",capsize=3)
    ax.set_xticks(range(len(labels)),labels,rotation=30);ax.set_ylabel("H-alpha pseudo-EW (Angstrom)")
    ax.set_title("One point per night; formal propagated noise only, not systematic uncertainty")
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/"night_EW.png",dpi=150);plt.close(fig)
    print(json.dumps({"primary":base,"beta_required":summary["continuum_only_beta_to_raise_2009_to_2008"],
                      "alpha_required":summary["broadband_alpha_to_dilute_2008_to_2009"],"MC":mc},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",required=True);p.add_argument("--out",required=True)
    a=p.parse_args();run(a.input,a.out)
