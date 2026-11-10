"""Exploratory CO-window shape sensitivity; not molecular identification."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from compare_spherex import source_filename
from compare_stellar_models import dust_grid


def residual_index(wave,observed,model):
    wave=np.asarray(wave);observed=np.asarray(observed);model=np.asarray(model)
    ref=(wave>=4.)&(wave<4.2);co=(wave>=4.55)&(wave<4.9)
    if ref.sum()<2 or co.sum()<2 or np.any(observed<=0) or np.any(model<=0):
        raise ValueError('Insufficient positive flux/window coverage')
    return float(np.median((observed/model)[co])/np.median((observed/model)[ref]))


def run(previous,models,dust,out):
    out.mkdir(parents=True,exist_ok=True)
    lookup={r['file']:r for r in json.loads((previous/'local/local_extraction.json').read_text())['records']}
    rows=[]
    for n in Table.read(previous/'native_IRSA_photometry.tbl',format='ascii.ipac'):
        r=lookup.get(source_filename(n['obs_publisher_did']))
        if r and r['primary_quality'] and int(n['flags'])==1<<21 and float(n['fit_ql'])<5:
            rows.append(dict(r,nwave=float(n['lambda']),nwidth=float(n['lambda_width']),nflux=float(n['flux']),nerror=float(n['flux_err'])))
    grain=dust_grid(dust/'draine_rv31_extinction.txt');results=[];retained=[]
    fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True,layout='constrained')
    for reduction,ax in zip(['local','native'],axes):
        for season in ['2025','2026']:
            a=[r for r in rows if r['season']==season]
            w=np.array([r['wavelength_um'] if reduction=='local' else r['nwave'] for r in a])
            width=np.array([r['width_um'] if reduction=='local' else r['nwidth'] for r in a])
            f=np.array([r['flux_ujy'] if reduction=='local' else r['nflux'] for r in a])
            x=(w[:,None]+width[:,None]*np.linspace(-.5,.5,31))*1e4
            opacity=np.exp(np.interp(np.log(x/1e4),np.log(grain[:,0]),np.log(grain[:,3]))-np.interp(np.log(.55),np.log(grain[:,0]),np.log(grain[:,3])))
            train=(w>=.75)&(w<=2.4);ref=(w>=4)&(w<4.2);co=(w>=4.55)&(w<4.9)
            best=None;accepted=[]
            for path in sorted(models.glob('ck*.fits')):
                with fits.open(path) as h:
                    wave=np.asarray(h[1].data['WAVELENGTH'],float)
                    for gravity in ['g10','g15','g20','g25','g30']:
                        fnu=np.asarray(h[1].data[gravity],float)*wave**2/2.99792458e18
                        base=np.interp(x,wave,fnu)
                        if np.any(base<=0):continue
                        for av in np.arange(0.,1.501,.05):
                            m=np.mean(base*10**(-.4*av*opacity),axis=1)
                            scale=float(np.exp(np.mean(np.log(f[train]/m[train]))));m*=scale
                            rms=float(np.sqrt(np.mean(np.log(f[train]/m[train])**2)))
                            case=dict(reduction=reduction,season=season,model=path.name,gravity=gravity,av=float(av),scale=scale,nearIR_log_rms=rms,CO_reference_residual_index=residual_index(w,f,m),
                                      reference_median_obs_model=float(np.median((f/m)[ref])),CO_median_obs_model=float(np.median((f/m)[co])))
                            if rms<=.05:accepted.append(case);retained.append(case)
                            if best is None or rms<best['nearIR_log_rms']:best=case;best_m=m.copy()
            summary=dict(reduction=reduction,season=season,n=len(a),n_ref=int(ref.sum()),n_co=int(co.sum()),best=best,compatible_cases=len(accepted),
                         compatible_CO_index_range=[min(r['CO_reference_residual_index'] for r in accepted),max(r['CO_reference_residual_index'] for r in accepted)] if accepted else None)
            results.append(summary)
            selection=w>=3.8;ratio=f/best_m/np.median((f/best_m)[ref])
            ax.plot(w[selection],ratio[selection],'.',label=season)
            ax.set_title(reduction+' reduction; best near-IR atmosphere, normalized at 4.0â€“4.2 Âµm')
        ax.axhline(1,c='gray');ax.axvspan(4.,4.2,alpha=.08,color='green');ax.axvspan(4.55,4.9,alpha=.08,color='orange')
        ax.set_ylabel('Observed/model, divided by reference');ax.legend();ax.grid(alpha=.2)
    axes[-1].set_xlabel('Wavelength (Âµm)');fig.savefig(out/'co_window_shape.png',dpi=160);plt.close(fig)
    report=dict(results=results,qualification='Post-hoc approximate rectangular responses; Rv3.1 only, coarse LTE plane-parallel CK04 logg1â€“3 grid. NearIR fit selects model; CO/reference is a broad window residual, not equivalent width or CO emission detection. 5% shape tolerance is nuisance sensitivity, not confidence interval. Native/local share pixels; detector/season correlated calibration and atmosphere error prevent assigning formal significance. No current spectral feature proves the cause of the 2018â€“2022 fading.')
    (out/'co_results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (out/'co_compatible_cases.json').write_text(json.dumps(retained,indent=2,allow_nan=False)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--previous',type=Path,required=True);p.add_argument('--models',type=Path,required=True);p.add_argument('--dust',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.previous,a.models,a.dust,a.out)
