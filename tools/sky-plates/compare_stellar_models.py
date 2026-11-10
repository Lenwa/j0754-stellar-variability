"""Retrospective g/W1 model fits and W2 checks; new-color predictions stay separate."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from fetch_stellar_grids import FILTERS
from stellar_models import read_bandpass,band_flux,magnitudes,interpolate,fit_temperature_scale,fit_partial_cover,fit_two_component,fit_spots,covering_predictions

BANDS=['g','W1','W2','B','V','r','I','J','H','Ks','i','z']
COLORS={'B-V':('B','V'),'g-r':('g','r'),'V-I':('V','I'),'J-H':('J','H'),'H-Ks':('H','Ks'),'g-Ks':('g','Ks')}


def dust_grid(path):
    rows=[]
    for line in path.read_text().splitlines():
        try:
            v=[float(x) for x in line.split()]
            if len(v)==6 and v[0]>0 and v[3]>0:rows.append(v)
        except ValueError:pass
    a=np.array(rows);return a[np.argsort(a[:,0])]


def extinction_deficits(wave,spectrum,bands,grain,tau):
    kref=np.interp(np.log(.48),np.log(grain[:,0]),np.log(grain[:,3]))
    columns=[]
    for b in bands:
        x=np.unique(np.r_[b.wave,wave[(wave>=b.wave.min())&(wave<=b.wave.max())]])
        w=np.interp(x,b.wave,b.transmission,left=0,right=0)*(x if b.detector else 1)
        flux=np.interp(x,wave,spectrum)
        opacity=np.exp(np.interp(np.log(x/1e4),np.log(grain[:,0]),np.log(grain[:,3]))-kref)
        transmitted=np.trapezoid(np.exp(-tau[:,None]*opacity)*flux*w,x,axis=1)/np.trapezoid(flux*w,x)
        columns.append(1-transmitted)
    return np.array(columns).T


def summarize(records):
    summary={}
    for family in sorted({r['family'] for r in records}):
        rows=[r for r in records if r['family']==family]
        good=[r for r in rows if r['max_training_residual']<=.02]
        def ranges(items,key):
            return {k:[float(min(r[key][k] for r in items)),float(max(r[key][k] for r in items))] for k in items[0][key]} if items else {}
        summary[family]=dict(cases=len(rows),training_compatible=len(good),
            best_training_residual=min(r['max_training_residual'] for r in rows),
            predicted_bands=ranges(good,'bands'),predicted_colors=ranges(good,'colors'),
            W2_residual_range=[min(r['W2_residual'] for r in good),max(r['W2_residual'] for r in good)] if good else None,
            qualification='Ranges of per-nuisance best fits passing g/W1 residual<=.02 only. For sampled full tolerance envelopes see compatibility_envelopes; neither is a confidence interval or W2-selected subset.')
    return summary


def main(inputs,dust,out):
    out.mkdir(exist_ok=True,parents=True)
    passbands={b:read_bandpass(inputs/name,fid) for b,(name,fid) in FILTERS.items()}
    temperatures=np.arange(3500.,6501.,250.)
    tau=np.geomspace(.001,300,800)
    grains={name:dust_grid(dust/f'draine_{name}_extinction.txt') for name in ['rv31','rv4','rv5.5']}
    records=[];coverage=[];envelopes={'primary':{},'sensitivity':{}}
    # These values come from previously saved matched halves, not a new fit to W2.
    seasons={'March':np.array([.312,.1175,.08875]),'October':np.array([.3225,.115,.09])}
    def save(family,result,meta,observed):
        if result is None:return
        pred=np.asarray(result['prediction']);result=dict(result);result.pop('prediction')
        records.append(dict(family=family,**meta,parameters=result,
            bands={b:float(p) for b,p in zip(BANDS,pred)},
            colors={c:float(pred[BANDS.index(a)]-pred[BANDS.index(b)]) for c,(a,b) in COLORS.items()},
            training_residuals=[float(x) for x in pred[:2]-observed[:2]],
            max_training_residual=float(np.max(np.abs(pred[:2]-observed[:2]))),
            W2_observed=float(observed[2]),W2_residual=float(pred[2]-observed[2])))
    def keep_envelope(family,predictions,meta,observed):
        if not len(predictions):return
        a=np.asarray(predictions)
        a=a[np.max(np.abs(a[:,:2]-observed[:2]),axis=1)<=.020000001]
        if not len(a):return
        branch='primary' if meta['response']=='g' and meta['season']=='March' else 'sensitivity'
        dest=envelopes[branch].setdefault(family,dict(sampled_predictions=0,bands={},colors={}))
        dest['sampled_predictions']+=len(a)
        for key,values in [('bands',{b:a[:,i] for i,b in enumerate(BANDS)}),
                           ('colors',{c:a[:,BANDS.index(x)]-a[:,BANDS.index(y)] for c,(x,y) in COLORS.items()})]:
            for name,v in values.items():
                bounds=[float(v.min()),float(v.max())]
                if name in dest[key]:bounds=[min(bounds[0],dest[key][name][0]),max(bounds[1],dest[key][name][1])]
                dest[key][name]=bounds
    for metal in ['ckm10','ckm05','ckp00']:
        native=[]
        for t in temperatures:
            with fits.open(inputs/f'{metal}_{int(t)}.fits') as h:
                wave=np.asarray(h[1].data['WAVELENGTH'],float)
                native.append({g:np.asarray(h[1].data[g],float) for g in ['g20','g25','g30']})
        for gravity in ['g20','g25','g30']:
            spectra=np.array([r[gravity] for r in native])
            for response in ['g','ztfg']:
                # ZTF response sensitivity is restricted to the nominal atmosphere.
                if response=='ztfg' and (metal!='ckm05' or gravity!='g25'):continue
                bands=[passbands[response] if b=='g' else passbands[b] for b in BANDS]
                grid=np.array([[band_flux(wave,s,b) for b in bands] for s in spectra])
                for bright_temp in [4000.,4500.,5000.,5500.]:
                    if response=='ztfg' and bright_temp!=5000:continue
                    bright=interpolate(temperatures,grid,bright_temp)
                    spec=interpolate(temperatures,spectra,bright_temp)
                    cool_t=np.arange(3500,bright_temp-249,10.)
                    cool_flux=np.array([interpolate(temperatures,grid,t) for t in cool_t])
                    area_t=np.arange(3500,bright_temp+1,2.)
                    area_flux=np.array([interpolate(temperatures,grid,t) for t in area_t])
                    cooling_delta=magnitudes(area_flux,bright)
                    screen={law:extinction_deficits(wave,spec,bands,grain,tau) for law,grain in grains.items()}
                    for season,observed in seasons.items():
                        meta=dict(metallicity=metal,logg=gravity,temperature_bright=bright_temp,response=response,season=season)
                        save('temperature_plus_scale',fit_temperature_scale(temperatures,grid,bright_temp,observed),meta,observed)
                        lo=np.max(observed[:2]-.02-cooling_delta[:,:2],axis=1)
                        hi=np.min(observed[:2]+.02-cooling_delta[:,:2],axis=1)
                        valid=lo<=hi
                        keep_envelope('temperature_plus_scale',np.r_[cooling_delta[valid]+lo[valid,None],cooling_delta[valid]+hi[valid,None]],meta,observed)
                        best=np.argmin(np.sum((cooling_delta[:,:2]-observed[:2])**2,axis=1))
                        save('fixed_area_cooling',dict(temperature_faint=float(area_t[best]),prediction=cooling_delta[best]),meta,observed)
                        keep_envelope('fixed_area_cooling',cooling_delta,meta,observed)
                        gray=float(np.mean(observed[:2]));save('gray_single_component',dict(prediction=np.full(len(BANDS),gray)),meta,observed)
                        for start in [0.,.1,.3]:
                            save('cool_spots',fit_spots(bright,cool_flux,cool_t,observed,start),meta,observed)
                            deficit=1-cool_flux/bright
                            allowed=covering_predictions(deficit/(1-start*deficit),observed,.02,1-start)
                            keep_envelope('cool_spots',[r['prediction'] for r in allowed],meta,observed)
                        for comp_temp in np.arange(3500,bright_temp,250):
                            comp=interpolate(temperatures,grid,comp_temp)
                            save('two_component_occultation',fit_two_component(bright,comp,observed),dict(meta,companion_temperature=float(comp_temp)),observed)
                            areas=np.geomspace(.001,100,240)
                            allowed=covering_predictions(1/(1+areas[:,None]*(comp/bright)),observed,.02)
                            keep_envelope('two_component_occultation',[r['prediction'] for r in allowed],meta,observed)
                        for law,deficit in screen.items():
                            save('partial_dust_covering',fit_partial_cover(tau,deficit,observed),dict(meta,dust_law=law),observed)
                            allowed=covering_predictions(deficit,observed,.02)
                            keep_envelope('partial_dust_covering',[r['prediction'] for r in allowed],meta,observed)
                            full=-2.5*np.log10(np.maximum(1-deficit,1e-250))
                            j=np.argmin(np.sum((full[:,:2]-observed[:2])**2,axis=1))
                            save('uniform_dust_screen',dict(tau=float(tau[j]),prediction=full[j]),dict(meta,dust_law=law),observed)
                            keep_envelope('uniform_dust_screen',full,meta,observed)
                    coverage.append(dict(metallicity=metal,logg=gravity,T=bright_temp,response=response))
    primary=[r for r in records if r['response']=='g' and r['season']=='March']
    sensitivity=[r for r in records if r not in primary]
    report=dict(protocol='protocol.json',training_bands=['g','W1'],check_band='W2',check_is_retrospective=True,
        observed={k:{b:float(v) for b,v in zip(BANDS,obs)} for k,obs in seasons.items()},
        summary=summarize(primary),sensitivity_summary=summarize(sensitivity),records=records,atmosphere_cases=coverage,
        compatibility_envelopes=envelopes,
        envelope_qualification='Sampled full g/W1 tolerance region (Â±.02mag each): tau800 logarithmic nodes; temperatures2K for cooling/scale and10K for spots; companion areas240lognodes .001â€“100. Analytic fraction/gray-offset interval endpoints included. No W2 selection; not a statistically calibrated or mathematically exhaustive confidence set.',
        limitations=['SDSSg approximates ASAS-SNg, and CK04 is a finite/coarse LTE spectral grid; atmospheric and response errors are not calibrated.',
                     'No statistical model probabilities or unique mechanism follows from fitting two training bands.',
                     'Nuisance ranges are not prior-weighted confidence intervals; boundary spot temperatures and alternative colder spots outside the grid matter.',
                     'Two-component fits are spectral dilution geometries, not validated binary orbits; geometric scale is not a measured radius change.',
                     'W2 was already known to the investigators but never used by fitting helpers; this is a retrospective held-out-band check.',
                     'Forecast colors refer to the fitted2018â†’2022 state change, not a prediction of the present2026 magnitude.'],
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [*[p for p in sorted(inputs.glob('*')) if p.suffix in ['.fits','.xml']],
                                *[dust/f'draine_{name}_extinction.txt' for name in grains],dust/'dust_analysis.json']})
    (out/'model_comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    families=['uniform_dust_screen','fixed_area_cooling','temperature_plus_scale','cool_spots','partial_dust_covering','two_component_occultation','gray_single_component']
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    for y,family in enumerate(families):
        s=report['summary'][family];env=envelopes['primary'].get(family)
        if env:
            axes[0].plot(env['bands']['W2'],[y,y],lw=5,solid_capstyle='round')
            axes[1].plot(env['colors']['V-I'],[y,y],lw=5,solid_capstyle='round')
        else:
            axes[0].text(.005,y,'Cannot fit g + W1 within 0.02 mag',va='center',fontsize=8)
    axes[0].axvline(.08875,c='black',ls='--',label='Measured W2 (excluded from fit)')
    axes[0].axvspan(.08875-.03,.08875+.03,color='gray',alpha=.13,label='Â±0.03 mag descriptive systematic allowance')
    labels=[x.replace('_',' ') for x in families]
    for ax in axes:ax.set_yticks(range(len(families)),labels);ax.set_ylim(len(families)-.5,-.5);ax.grid(axis='x',alpha=.2)
    axes[0].set_xlim(-.025,.21);axes[0].set_xlabel('Predicted W2 fading (mag)')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=2,fontsize=8)
    axes[1].set_xlabel('Predicted increase in V âˆ’ I (mag)');axes[1].set_title('Additional-color predictions, currently untested')
    fig.suptitle('J0754: competing explanations constrained by g and W1 only\nSampled atmosphere/geometry regions within Â±0.02 mag in both training bands; not confidence intervals')
    fig.savefig(out/'model_predictions.png',dpi=170);plt.close(fig)
    for family,s in report['summary'].items():
        env=envelopes['primary'].get(family,{})
        print(family,s['training_compatible'],'/',s['cases'],'best-fit W2',s['predicted_bands'].get('W2'),
              'sampled compatible W2',env.get('bands',{}).get('W2'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--dust',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.inputs,a.dust,a.out)
