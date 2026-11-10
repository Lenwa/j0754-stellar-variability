"""Exploratory, within-grating Halpha blue-wing check on preserved FAST FITS."""
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar, curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fast_line_tests import load_spectra

OUT=Path('results/sky-plates-j0754-fast-bluewing-20261108')
INPUT=Path('results/sky-plates-j0754-spherex-20261108/literature/fast')


def normalize(r):
    x,y=r['x'],r['y']
    m=((x>=6520)&(x<=6540))|((x>=6585)&(x<=6610))
    c=np.polyval(np.polyfit(x[m]-6562.8,y[m],1),x-6562.8)
    return y/c,np.sqrt(r['var'])/c


def main():
    if not (OUT/'protocol.txt').exists():raise RuntimeError('Preserve exploratory protocol first')
    rows=load_spectra(INPUT)
    for r in rows:r['f'],r['err']=normalize(r)
    fig,axs=plt.subplots(2,2,figsize=(12,8))
    results=[]
    for group,g in enumerate(['300','600']):
        rr=[r for r in rows if r['grating']==g];ref=rr[0]
        x=ref['x'];metal=((x>=6400)&(x<=6500))|((x>=6600)&(x<=6720))
        # High-pass the counts spectrum to remove continuum colour/slit throughput.
        for r in rr:
            r['hp']=r['y']/gaussian_filter1d(r['y'],20/np.median(np.diff(r['x']))) - 1
        def loss(shift,r):
            return np.mean((np.interp(x[metal]+shift,r['x'],r['hp'])-ref['hp'][metal])**2)
        aligned=[]
        grid=np.arange(6538,6588,.25)
        for r in rr:
            shift=minimize_scalar(lambda s:loss(s,r),bounds=(-3,3),method='bounded').x
            r['shift']=float(shift)
            aligned.append(np.interp(grid+shift,r['x'],r['f']))
        # Reference core defines relative velocity zero; not an absolute RV.
        m=(grid>6557)&(grid<6569)
        def core(x,b,a,mu,s):return b-a*np.exp(-.5*((x-mu)/s)**2)
        fit,_=curve_fit(core,grid[m],aligned[0][m],p0=[1,.15,6563,2],
                        bounds=([.8,0,6560,.5],[1.2,.6,6566,7]))
        center=float(fit[2]);vel=(grid-center)/center*299792.458
        # Equal night weighting: three repeats on one night are not three epochs.
        nights=sorted({r['night'] for r in rr})
        means={n:np.mean([f for r,f in zip(rr,aligned) if r['night']==n],axis=0) for n in nights}
        for n in nights:
            z=means[n];others=np.mean([means[k] for k in nights if k!=n],axis=0)
            # Fit smooth continuum and symmetric core-depth changes, excluding
            # the candidate blue interval, to prevent broad fill-in being called gas.
            sym=(others+np.interp(2*center-grid,grid,others))/2
            design=np.c_[np.ones(len(grid)),(grid-center)/20,1-sym]
            blue=(vel>=-600)&(vel<=-150)
            mask=(np.abs(vel)<900)&~blue
            beta=np.linalg.lstsq(design[mask],(z-others)[mask],rcond=None)[0]
            residual=z-others-design@beta
            ew=float(-np.trapezoid(residual[blue],grid[blue]))
            peak=float(-min(residual[blue]))
            results.append(dict(night=n,grating=g,reference_core_A=center,
                shifts_A=[r['shift'] for r in rr if r['night']==n],
                blue_interval_km_s=[-600,-150],extra_absorption_EW_A=ew,
                maximum_blue_residual_depth_fraction=peak,
                note='Descriptive residual; continuum, interpolation, template and LSF systematics not converted into discovery significance.'))
            axs[group,0].plot(vel,z,label=n,lw=1)
            axs[group,1].plot(vel,residual,label=n,lw=1)
        for col in [0,1]:
            ax=axs[group,col];ax.axvspan(-600,-150,color='skyblue',alpha=.15)
            ax.axvline(0,color='.5',lw=.6);ax.set_xlim(-700,700);ax.legend(fontsize=7)
            ax.set_xlabel('Velocity relative to reference Halpha core (km/s)')
            ax.grid(alpha=.12)
        axs[group,0].set(title=f'FAST {g}: registered night spectra',ylabel='Flux / local continuum')
        axs[group,1].axhline(0,color='black',lw=.7)
        axs[group,1].set(title='Residual after other-night template + core-depth change',ylabel='Extra flux / continuum')
        axs[group,1].set_ylim(-.08,.08)
    fig.suptitle('Search for blue Halpha absorption â€” shaded: -600 to -150 km/s')
    fig.tight_layout();fig.savefig(OUT/'registered_Halpha_bluewing.png',dpi=160)
    result=dict(input=str(INPUT),spectra=len(rows),independent_nights=len({r['night'] for r in rows}),
        assumed_native_FWHM_A={'300':7,'600':3},
        corresponding_resolution_km_s={'300':7/6562.8*299792.458,'600':3/6562.8*299792.458},
        method='Native spectra, same-grating comparisons, neighbouring metal-line registration, equal-night reference; smooth continuum and symmetric line-depth nuisance. Exposures retained separately in input; plotted/compared night means.',
        assessment='No clean detached or consistently identifiable blue absorption component in native or registered spectra. Low-resolution 2008-10-22 has a blue/core residual near -200km/s, blended with the main line at assumed ~320km/s resolution; not classified as an eruption. Higher-resolution blue-side dips recur at the same relative wavelengths and do not show a distinct evolving component.',
        limits='No precise absolute RV or calibrated LSF. No atmosphere model or injected detection limits; fixed blends are not eruptions. Sparse 2008-2012 spectra are not contemporaneous with 2021 TESS flare. Weak/blended or unobserved eruptions remain unconstrained.',records=results)
    (OUT/'bluewing_results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
