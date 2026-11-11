"""Offline photometric tests of Tang's other two historical fading giants."""
import hashlib
import json
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis_core import separation_arcsec
from screen_giant_analogs import ROOT,write_csv,gap_series
from tang_sibling_core import absolute_time,quality_mask,bin_series,periodic_fit,scan_periods,block_sign_diagnostic,excursion_runs,local_residual_scale,clean

INPUT=ROOT/'results/giant-screening-20261111'
OUT=ROOT/'results/tang-sibling-analysis-20261111'
TARGETS=dict(J0736=dict(period=39.42091782293346,ra=114.02711824252,dec=21.23634538304,tic=70330439),J0830=dict(period=38.38919173490344,ra=127.66033651549,dec=14.12030855354,tic=20445282,epic=211569759))
PERIODS=1/np.linspace(1/100,1/5,1601)


def read_inputs(write_outputs=True):
    records=[];curves=defaultdict(list);native={}
    for product in json.loads((INPUT/'sibling_lightcurve_products.json').read_text()):
        path=INPUT/product['file'];target=product['target'];expected=TARGETS[target]
        with fits.open(path) as h:
            d=h[1].data;header=h[1].header;primary=h[0].header;mission='K2' if 'KEPLERID' in primary else 'TESS'
            identity=int(primary['KEPLERID'] if mission=='K2' else primary['TICID'])
            if identity!=expected['epic' if mission=='K2' else 'tic']:raise ValueError('Wrong stellar identity')
            separation=separation_arcsec(primary['RA_OBJ'],primary['DEC_OBJ'],expected['ra'],expected['dec'])
            if separation>1:raise ValueError('Header position exceeds one arcsecond match')
            t=absolute_time(d['TIME'],header);q=np.asarray(d['SAP_QUALITY' if mission=='K2' else 'QUALITY'])
            group='K2C5' if mission=='K2' else 'S'+str(primary['SECTOR'])
            columns=['SAP_FLUX','PDCSAP_FLUX'] if mission=='K2' else ['SAP_FLUX','KSPSAP_FLUX','KSPSAP_FLUX_SML','KSPSAP_FLUX_LAG']
            extra_names=['POS_CORR1','POS_CORR2','SAP_BKG'] if mission=='K2' else ['SAP_X','SAP_Y','SAP_BKG']
            raw_extra=np.column_stack([np.asarray(d[name],float) for name in extra_names])
            native[path.name]=dict(target=target,mission=mission,t=t,quality=q,cadence=np.asarray(d['CADENCENO']),data={name:np.asarray(d[name],float) for name in columns},extra=raw_extra,group=group)
            for quality in ['strict']+(['kepler_default'] if mission=='K2' else []):
                for col in columns:
                    flux=np.asarray(d[col],float);mask=quality_mask(q,quality)&np.isfinite(t)&np.isfinite(flux)&(flux>0)
                    tt,ff=t[mask],flux[mask];center=float(np.median(ff));y=ff/center-1
                    extra=raw_extra[mask].copy();imputed=[]
                    for j in range(extra.shape[1]):
                        valid=np.isfinite(extra[:,j]);imputed.append(int(sum(~valid)))
                        extra[~valid,j]=float(np.median(extra[valid,j])) if sum(valid) else 0.
                    g=np.full(len(tt),group);binned=bin_series(tt,y,g,extra=extra)
                    key=(target,mission,col,quality);curves[key].append(binned)
                    correlations={name:float(np.corrcoef(y,extra[:,j])[0,1]) if np.std(extra[:,j])>0 else None for j,name in enumerate(extra_names)}
                    record=dict(file=path.name,target=target,mission=mission,column=col,quality=quality,group=group,raw_rows=len(d),accepted_rows=len(tt),binned_rows=len(binned['t']),median_native_flux=center,first_bjd=float(min(tt)),last_bjd=float(max(tt)),span_days=float(np.ptp(tt)),p05_p95_relative_percent=((np.quantile(ff/center,[.05,.95])-1)*100).tolist(),header_match_arcsec=separation,header_identifier=identity,regressor_imputed_rows=imputed,correlations=correlations,camera=primary.get('CAMERA'),ccd=primary.get('CCD'))
                    records.append(record)
                    if write_outputs:write_csv(OUT/f'{target}_{mission}_{group}_{col}_{quality}_bins.csv',[dict(bjd=float(ti),relative_flux=float(yi),group=str(gi),count=int(ci),x=float(e[0]),y_position=float(e[1]),background=float(e[2])) for ti,yi,gi,ci,e in zip(binned['t'],binned['y'],binned['group'],binned['count'],binned['extra'])])
    combined={key:{name:np.concatenate([item[name] for item in values]) for name in ['t','y','group','extra','count']} for key,values in curves.items()}
    return records,combined,native


def common_cadences(native):
    result=[]
    for group in ['S44','S45','S46']:
        a=next(r for r in native.values() if r['target']=='J0736' and r['group']==group)
        b=next(r for r in native.values() if r['target']=='J0830' and r['group']==group)
        for col in ['SAP_FLUX','KSPSAP_FLUX']:
            ma=(a['quality']==0)&np.isfinite(a['t'])&np.isfinite(a['data'][col]);mb=(b['quality']==0)&np.isfinite(b['t'])&np.isfinite(b['data'][col])
            _,ia,ib=np.intersect1d(a['cadence'][ma],b['cadence'][mb],return_indices=True)
            t=(a['t'][ma][ia]+b['t'][mb][ib])/2
            ya=a['data'][col][ma][ia];yb=b['data'][col][mb][ib];ya=ya/np.median(ya)-1;yb=yb/np.median(yb)-1
            bins=bin_series(t,ya,np.zeros(len(t)),extra=yb[:,None]);x=bins['t']-np.mean(bins['t'])
            detrended=[]
            for values in [bins['y'],bins['extra'][:,0]]:detrended.append(values-np.polyval(np.polyfit(x,values,1),x))
            result.append(dict(group=group,column=col,matched_cadences=len(t),binned_pairs=len(x),linear_detrended_correlation=float(np.corrcoef(*detrended)[0,1]),interpretation='Different variable stars on common spacecraft cadences; correlation is a systematics diagnostic, not a constant-star calibration.'))
    return result


def preliminary_excursions(native):
    records=[]
    for file,r in native.items():
        col='PDCSAP_FLUX' if r['mission']=='K2' else 'KSPSAP_FLUX'
        flux=r['data'][col];good=(r['quality']==0)&np.isfinite(r['t'])&np.isfinite(flux)&(flux>0)
        t=r['t'][good];y=flux[good]/np.median(flux[good])-1
        cadence_numbers=r['cadence'][good]
        order=np.argsort(t);t,y,cadence_numbers=t[order],y[order],cadence_numbers[order]
        baseline=np.array([np.median(y[np.abs(t-ti)<=.5]) for ti in t]);residual=y-baseline
        sigma=local_residual_scale(t,residual)
        cadence=float(np.median(np.diff(r['t'][np.isfinite(r['t'])])))
        cadence_step=int(np.median(np.diff(r['cadence'])))
        runs=excursion_runs(t,residual>5*sigma,cadence,cadence_numbers=cadence_numbers,cadence_step=cadence_step)
        for run in runs:
            peak=run[int(np.argmax(residual[run]))];amplitude=float(residual[peak]);limit=.005 if r['mission']=='K2' else .002
            if amplitude<limit:continue
            records.append(dict(file=file,target=r['target'],mission=r['mission'],group=r['group'],column=col,peak_bjd=float(t[peak]),utc_date=Time(t[peak],format='jd',scale='tdb').utc.isot,peak_relative_flux=amplitude,local_residual_sigma=float(sigma[peak]),noise_half_window_days=.5,noise_minimum_points=20,consecutive_points=len(run),native_cadence_step=cadence_step,accepted_span_hours=float((t[run[-1]]-t[run[0]])*24),status='Unvetted positive excursion; no pixel localization or stellar-flare classification',_t=t,_y=y,_baseline=baseline))
    return records


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    protocol=OUT/'protocol.json';plan=ROOT/'docs/superpowers/plans/2026-11-11-tang-siblings.md'
    if not protocol.exists():
        snapshot=OUT/'protocol_plan_snapshot.txt';snapshot.write_bytes(plan.read_bytes())
        protocol.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),status='Exploratory follow-up; native headers and flux ranges inspected, period-fit results not yet seen.',plan=plan.relative_to(ROOT).as_posix(),plan_original_snapshot=snapshot.relative_to(ROOT).as_posix(),plan_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),period_grid_days=[5,100],frequency_grid_points=len(PERIODS),quality_primary='All flags zero; K2 default mask as separate sensitivity',period_hypotheses={k:v['period'] for k,v in TARGETS.items()},independent_epochs='K2 2015 and TESS 2021 fitted independently; alternative reductions not independent instruments',documentation=['https://archive.stsci.edu/hlsp/qlp','https://lightkurve.github.io/lightkurve/_modules/lightkurve/io/qlp.html','https://lightkurve.github.io/lightkurve/_modules/lightkurve/utils.html']),indent=2)+'\n',encoding='utf8')
    records,curves,native=read_inputs();results=[];fits_for_plot={};scans_for_plot={}
    for key,data in curves.items():
        target,mission,col,quality=key;period=TARGETS[target]['period'];t,y,g,extra=[data[n] for n in ['t','y','group','extra']]
        label='_'.join(key);entry=dict(target=target,mission=mission,column=col,quality=quality,points=len(t),span_days=float(np.ptp(t)),models={})
        modes=['global_linear','sector_linear','sector_quadratic','instrument'] if col in ['SAP_FLUX','KSPSAP_FLUX','PDCSAP_FLUX'] else ['global_linear']
        for mode in modes:
            one=periodic_fit(t,y,g,period,mode,extra=extra);two=periodic_fit(t,y,g,period,mode,harmonics=2,extra=extra)
            item=dict(at_gaia=clean(one),at_half_gaia=clean(periodic_fit(t,y,g,period/2,mode,extra=extra)),at_gaia_two_harmonics=clean(two))
            if col in ['SAP_FLUX','KSPSAP_FLUX','PDCSAP_FLUX']:
                peak,power=scan_periods(t,y,g,PERIODS,mode,extra=extra);item['free_peak']=peak
                write_csv(OUT/f'{label}_{mode}_periodogram.csv',[dict(period_days=float(p),power=float(v)) for p,v in zip(PERIODS,power)])
                if mode in ['global_linear','sector_linear'] and quality=='strict':
                    scans_for_plot[(target,mission,col,mode)]=power
                    if col=='SAP_FLUX':item['block_sign_diagnostics']=block_sign_diagnostic(t,y,g,PERIODS,mode,extra=extra)
                if mode=='global_linear':
                    peak2,power2=scan_periods(t,y,g,PERIODS,mode,harmonics=2,extra=extra);item['free_two_harmonic_peak']=peak2
                    write_csv(OUT/f'{label}_{mode}_two_harmonic_periodogram.csv',[dict(period_days=float(p),power=float(v)) for p,v in zip(PERIODS,power2)])
                fits_for_plot[(target,mission,col,quality,mode)]=one
            entry['models'][mode]=item
        results.append(entry)
        print(label,'Gaia power',round(entry['models']['global_linear']['at_gaia']['power'],3),'peak',entry['models']['global_linear'].get('free_peak',{}).get('best_period_days'),flush=True)
    # Equal-epoch mean power permits independent phases/amplitudes. The baseline
    # remains each ~75-day epoch, never the six-year distance between them.
    combined=[]
    for mode in ['global_linear','sector_linear']:
        a=scans_for_plot[('J0830','TESS','SAP_FLUX',mode)];b=scans_for_plot[('J0830','K2','SAP_FLUX',mode)];v=(a+b)/2;i=int(np.argmax(v))
        combined.append(dict(target='J0830',mode=mode,best_period_days=float(PERIODS[i]),mean_epoch_power=float(v[i]),at_gaia=float(np.interp(1/TARGETS['J0830']['period'],1/PERIODS,v)),phase_coherence='Not imposed across epochs',interpretation='Mean descriptive profile, not a six-year coherent period measurement'))
    shared=common_cadences(native);excursions=preliminary_excursions(native)
    result=dict(file_measurements=records,period_tests=results,independent_epoch_profiles=combined,common_cadence_diagnostics=shared,preliminary_excursions=clean(excursions),limitations=['Gaia SB1 fits remain candidates.','Only about two proposed cycles per observing epoch.','QLP detrending removes slow structure; corrected curves cannot set a stellar-amplitude limit.','SAP is less filtered but retains instrumental effects.','Sector-normalized data cannot measure absolute decade-long fading.','No pixel-level flare host attribution in this stage.'])
    (OUT/'results.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n',encoding='utf8')
    write_csv(OUT/'file_measurements.csv',records)
    plots(curves,fits_for_plot,scans_for_plot,shared,excursions)
    hashes={str((INPUT/name).relative_to(ROOT).as_posix()):hashlib.sha256((INPUT/name).read_bytes()).hexdigest() for name in native}
    for source in [Path(__file__),Path(__file__).with_name('tang_sibling_core.py')]:hashes[source.relative_to(ROOT).as_posix()]=hashlib.sha256(source.read_bytes()).hexdigest()
    (OUT/'input_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf8')
    print('EXCURSIONS',len(excursions),flush=True)


def plots(curves,fitmap,scanmap,shared,excursions):
    epochs=[('J0736','TESS'),('J0830','TESS'),('J0830','K2')]
    fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
    for i,(target,mission) in enumerate(epochs):
        raw=curves[(target,mission,'SAP_FLUX','strict')];col='KSPSAP_FLUX' if mission=='TESS' else 'PDCSAP_FLUX';corrected=curves[(target,mission,col,'strict')]
        for ax in axes[i]:ax.grid(alpha=.15)
        for group in sorted(set(raw['group'])):
            mask=raw['group']==group;t=raw['t'][mask];x=t-min(raw['t']);axes[i,0].plot(*gap_series(x,raw['y'][mask]*100,.5),'.-',ms=2,lw=.5,label=group)
            fit=fitmap[(target,mission,'SAP_FLUX','strict','global_linear')]
            axes[i,0].plot(*gap_series(x,(fit['_periodic'][mask]+fit['_nuisance'][mask])*100,.5),color='black',lw=1)
        axes[i,0].set(title=target+' / '+mission+' aperture flux and fixed-Gaia fit',xlabel='Days within this epoch',ylabel='Relative flux (%)');axes[i,0].legend(fontsize=7)
        for group in sorted(set(corrected['group'])):
            m=corrected['group']==group;axes[i,1].plot(*gap_series(corrected['t'][m]-min(corrected['t']),corrected['y'][m]*100,.5),'.-',ms=2,lw=.5,label=group)
        axes[i,1].set(title=target+' / '+mission+' pipeline-corrected flux',xlabel='Days within this epoch',ylabel='Relative flux (%)')
    fig.suptitle('Separate epochs and reductions: fixed Gaia-period fits are descriptive\nIndependent sector medians; slow filtering and instrumental structure remain visible.')
    fig.savefig(OUT/'sibling_lightcurves.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,10),layout='constrained')
    for ax,(target,mission) in zip(axes,epochs):
        col='KSPSAP_FLUX' if mission=='TESS' else 'PDCSAP_FLUX'
        for mode,color in [('global_linear','navy'),('sector_linear','darkorange')]:
            ax.plot(PERIODS,scanmap[(target,mission,'SAP_FLUX',mode)],color=color,label='SAP / '+mode)
            ax.plot(PERIODS,scanmap[(target,mission,col,mode)],color=color,ls=':',alpha=.7,label=col+' / '+mode)
        p=TARGETS[target]['period'];ax.axvline(p,color='black',ls='--',label=f'Gaia {p:.2f} d');ax.axvline(p/2,color='gray',ls='--',label='Half Gaia period')
        ax.set(title=target+' / '+mission,xlim=(5,100),ylim=(0,1),ylabel='Fractional variance improvement',xlabel='Trial period (days)');ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.15)
    fig.suptitle('Period profiles depend on baseline and reduction; peaks are not significance levels')
    fig.savefig(OUT/'period_sensitivity.png',dpi=160);plt.close(fig)
    if excursions:
        fig,axes=plt.subplots(len(excursions),1,figsize=(9,max(3,2.5*len(excursions))),layout='constrained',squeeze=False)
        for ax,event in zip(axes[:,0],excursions):
            t=event['_t'];m=np.abs(t-event['peak_bjd'])<1;ax.plot((t[m]-event['peak_bjd'])*24,event['_y'][m]*100,'.',ms=2);ax.plot((t[m]-event['peak_bjd'])*24,event['_baseline'][m]*100,lw=.8)
            ax.set(title=event['target']+' '+event['utc_date']+' / unvetted excursion',xlabel='Hours from peak',ylabel='Corrected relative flux (%)');ax.grid(alpha=.2)
        fig.savefig(OUT/'unvetted_excursions.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
