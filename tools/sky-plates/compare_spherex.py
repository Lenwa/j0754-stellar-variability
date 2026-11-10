"""Exploratory static SED check, retaining reduction and calibration sensitivity.

Rectangular CWAVE/CBAND responses are approximations, not the native passbands.
No model probabilities, confidence limits or state-change causal claims.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from compare_stellar_models import dust_grid


def source_filename(identifier):
    return str(identifier).split('?')[-1].replace('/','')+'.fits'


def rectangular_fnu(wave,flambda,centers,widths):
    x=(np.asarray(centers)[:,None]+np.asarray(widths)[:,None]*np.linspace(-.5,.5,31))*1e4
    # Interpolate F_nu rather than F_lambda, preserving constant F_nu exactly.
    fnu=np.asarray(flambda)*np.asarray(wave)**2/2.99792458e18
    return np.mean(np.interp(x,wave,fnu),axis=1)


def main(root,models,dust):
    local=json.loads((root/'local/local_extraction.json').read_text())
    native=Table.read(root/'native_IRSA_photometry.tbl',format='ascii.ipac')
    lookup={r['file']:r for r in local['records']}
    rows=[]
    for a in native:
        r=lookup.get(source_filename(a['obs_publisher_did']))
        if r is None:continue
        rows.append(dict(r,native_flux=float(a['flux']),native_error=float(a['flux_err']),
                         native_flags=int(a['flags']),native_wave=float(a['lambda']),
                         native_width=float(a['lambda_width']),native_fit_ql=float(a['fit_ql'])))
    good=[r for r in rows if r['primary_quality']]
    # Additional official flag restriction is a sensitivity, not retrospective
    # replacement of the predeclared local quality rule. SOURCE is informational.
    clean=[r for r in good if r['native_flags']==1<<21 and r['native_fit_ql']<5]
    comparison={}
    for d in range(1,7):
        a=[r for r in good if r['detector']==d]
        comparison[str(d)]={'n':len(a),'native_local_flux_ratio_quantiles':np.quantile([r['native_flux']/r['flux_ujy'] for r in a],[0,.05,.5,.95,1]).tolist(),
                            'native_minus_local_wavelength_um':np.quantile([r['native_wave']-r['wavelength_um'] for r in a],[0,.5,1]).tolist(),
                            'median_reduced_chi2':float(np.median([r['reduced_chi2'] for r in a]))}
    grains={rv:dust_grid(dust/f'draine_{rv}_extinction.txt') for rv in ['rv31','rv4','rv5.5']}
    cases=[];best_curves={}
    for reduction in ['local','native']:
        for season in ['2025','2026']:
            a=[r for r in clean if r['season']==season]
            center=np.array([r['wavelength_um'] if reduction=='local' else r['native_wave'] for r in a])
            width=np.array([r['width_um'] if reduction=='local' else r['native_width'] for r in a])
            observed=np.array([r['flux_ujy'] if reduction=='local' else r['native_flux'] for r in a])
            train=(center>=.75)&(center<=2.4)
            x=(center[:,None]+width[:,None]*np.linspace(-.5,.5,31))*1e4
            best=None;compatible=[]
            for path in sorted(models.glob('ck*.fits')):
                with fits.open(path) as h:
                    wave=np.asarray(h[1].data['WAVELENGTH'],float)
                    for gravity in ['g20','g25','g30']:
                        fnu=np.asarray(h[1].data[gravity],float)*wave**2/2.99792458e18
                        spectrum=np.interp(x,wave,fnu)
                        if np.any(spectrum<=0):continue
                        for rv,grain in grains.items():
                            opacity=np.exp(np.interp(np.log(x/1e4),np.log(grain[:,0]),np.log(grain[:,3]))-np.interp(np.log(.55),np.log(grain[:,0]),np.log(grain[:,3])))
                            for av in np.arange(0.,2.001,.05):
                                prediction=np.mean(spectrum*10**(-.4*av*opacity),axis=1)
                                scale=np.exp(np.mean(np.log(observed[train]/prediction[train])))
                                prediction*=scale
                                score=float(np.mean(np.log(observed[train]/prediction[train])**2))
                                if score<=.05**2:
                                    compatible.append([float(np.median((observed/prediction)[(center>=lo)&(center<hi)]))
                                                       for lo,hi in [(2.45,3.75),(3.8,4.4),(4.4,5.)]])
                                if best is None or score<best['log_rms']**2:
                                    best=dict(model=path.name,gravity=gravity,rv=rv,av=float(av),scale=float(scale),log_rms=float(np.sqrt(score)))
                                    best_prediction=prediction.copy()
            ratios=observed/best_prediction
            best['five_percent_log_rms_sensitivity']={'sampled_cases':len(compatible),
                'MIR_median_ratio_bounds':np.column_stack((np.min(compatible,axis=0),np.max(compatible,axis=0))).tolist() if compatible else [],
                'qualification':'Post-hoc loose shape tolerance, not confidence bounds; includes temperature/gravity/metallicity/extinction degeneracy in sampled grid.'}
            for lo,hi in [(2.45,3.75),(3.8,4.4),(4.4,5.0)]:
                selection=(center>=lo)&(center<hi)
                best[f'ratio_{lo}_{hi}']=dict(n=int(selection.sum()),median=float(np.median(ratios[selection])),range=np.quantile(ratios[selection],[.05,.95]).tolist())
            cases.append(dict(reduction=reduction,season=season,n=int(train.sum()),**best))
            best_curves[(reduction,season)]=(center,observed,best_prediction)
    report=dict(quality_rows=len(good),joint_clean_rows=len(clean),detector_comparison=comparison,photosphere_cases=cases,
                qualification='Post-hoc exploratory shape check. Fits only .75â€“2.4um in each recovery-era season; MIR already seen before model calculation, so retrospective prediction check. Rectangular CWAVE/CBAND responses, coarse CK04 grid (3500â€“6500K/logg2â€“3/metalâˆ’1â€“0), Draine static foreground attenuation and free normalization. Equal log-flux weights; grid score is not likelihood. Parameter degeneracy, correlated calibration, missing actual response curves and atmosphere errors prevent precision temperatures or excess limits. Local and native methods share science pixels; their agreement is methodological, not independent astrophysical replication.')
    (root/'spectrum_comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True,layout='constrained')
    for (reduction,season),(center,observed,prediction) in best_curves.items():
        if reduction!='local':continue
        order=np.argsort(center)
        axes[0].plot(center,observed/1000,'.',label=season+' measured')
        axes[0].plot(center[order],prediction[order]/1000,'-',alpha=.7,label=season+' photosphere fitted below 2.4um')
        axes[1].plot(center,observed/prediction,'.',label=season)
    axes[0].set_ylabel('Flux density (mJy)');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Measured / photosphere');axes[1].axhline(1,c='gray');axes[1].axvline(2.4,c='gray',ls=':');axes[1].legend()
    axes[1].set_xlabel('Wavelength (Âµm)')
    axes[0].set_title('Exploratory stellar-atmosphere check; approximate rectangular responses')
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(root/'photosphere_check.png',dpi=160);plt.close(fig)
    print(json.dumps(cases,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--models',type=Path,required=True);p.add_argument('--dust',type=Path,required=True)
    a=p.parse_args();main(a.root,a.models,a.dust)
