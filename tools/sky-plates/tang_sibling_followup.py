"""Explicit post-initial-fit robustness tests; no independent observations implied."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_tang_siblings import ROOT,OUT,INPUT,TARGETS,PERIODS,read_inputs
from tang_sibling_core import periodic_fit,heldout_score,scan_periods,block_sign_diagnostic,clean,absolute_time,quality_mask,bin_series,restore_k2sc_variability,arc_labels
from screen_giant_analogs import write_csv,gap_series


def fetch_alternatives():
    import fetch_giant_screening as retrieval
    retrieval.OUT=OUT
    metadata=json.loads((INPUT/'J0830_mast.json').read_text())
    selected=[]
    for row in metadata['data']:
        uri=str(row.get('dataURL') or '')
        if not uri.startswith(('mast:HLSP/everest/','mast:HLSP/k2sc/')):continue
        file=uri.rsplit('/',1)[-1];receipt_path=OUT/(file+'.receipt.json')
        if receipt_path.exists() and not (OUT/file).exists():receipt=json.loads(receipt_path.read_text())
        else:receipt=retrieval.fetch(file,'https://mast.stsci.edu/api/v0.1/Download/file',dict(uri=uri))
        selected.append(dict(target='J0830',file=file,uri=uri,state=receipt['state'],reason='Post-initial-fit comparison; alternate correction of the same K2 epoch, not an independent observation'))
    (OUT/'additional_k2_products.json').write_text(json.dumps(selected,indent=2)+'\n',encoding='utf8')


def robustness(curves):
    outputs=[];corrected=[]
    for target in TARGETS:
        data=curves[(target,'TESS','SAP_FLUX','strict')];t,y,g,extra=[data[k] for k in ['t','y','group','extra']];p=TARGETS[target]['period']
        record=dict(target=target,followup_status='Specified after initial results showed a 38.4-day instrument-regression peak in J0830',regressor_subsets=[],leave_one_sector_out=[])
        for name,indices in [('x_only',[0]),('y_only',[1]),('background_only',[2]),('xy',[0,1]),('x_background',[0,2]),('y_background',[1,2]),('all',[0,1,2])]:
            z=extra[:,indices];fit=periodic_fit(t,y,g,p,'instrument',extra=z);peak,power=scan_periods(t,y,g,PERIODS,'instrument',extra=z)
            record['regressor_subsets'].append(dict(regressors=name,at_gaia=clean(fit),free_peak=peak))
        record['conditional_block_sign']=block_sign_diagnostic(t,y,g,PERIODS,'instrument',extra=extra)
        full=periodic_fit(t,y,g,p,'instrument',extra=extra)
        for sector in sorted(set(g)):
            train=g!=sector;test=g==sector
            fit=periodic_fit(t[train],y[train],g[train],p,'instrument',extra=extra[train])
            held=heldout_score(t[test],y[test],g[test],p,'instrument',fit['beta'],extra=extra[test])
            record['leave_one_sector_out'].append(dict(sector=sector,training_beta=fit['beta'],training_amplitude_percent=fit['harmonic_semiamplitudes'][0]*100,heldout=held))
        # Corrected curve retains the fitted periodic part; instrument coefficients
        # were estimated jointly, not from a pre-subtracted time series.
        corrected.append(dict(target=target,t=t,y=y-full['_nuisance'],g=g,fit=full['_periodic'],period=p))
        outputs.append(record)
    return outputs,corrected


def k2sc():
    path=OUT/'hlsp_k2sc_k2_llc_211569759-c05_kepler_v2_lc.fits';records=[];plots=[]
    if not path.exists():return records,plots
    with fits.open(path) as h:
        for extension in ['PDC','SAP']:
            d=h[extension].data;header=h[extension].header
            if int(header['EPIC'])!=TARGETS['J0830']['epic']:raise ValueError('Wrong K2SC EPIC identity')
            t=absolute_time(d['time'],header);flux=restore_k2sc_variability(d['flux'],d['trtime'])
            for quality in ['strict','kepler_default']:
                # Permit the redundant K2-quality marker only in the relaxed case;
                # reject K2SC outlier/nonfinite flags in both cases.
                mflag=np.asarray(d['mflags'],np.int64);allowed=(mflag==0) if quality=='strict' else ((mflag&~1)==0)
                mask=quality_mask(d['quality'],quality)&allowed&np.isfinite(t)&np.isfinite(flux)&(flux>0)
                tt=t[mask];yy=flux[mask]/np.median(flux[mask])-1
                b=bin_series(tt,yy,np.zeros(len(tt)));g=b['group'];p=TARGETS['J0830']['period']
                record=dict(extension=extension,quality=quality,native_accepted=len(tt),binned_points=len(b['t']),span_days=float(np.ptp(tt)),upstream_data_release=int(h[0].header['DATA_REL']),formula='flux + trtime - median(trtime); restores the time component as recommended by K2SC',time_component_column='trtime',position_component_column='trposi',models={})
                for mode in ['global_linear','sector_quadratic']:
                    fit=periodic_fit(b['t'],b['y'],g,p,mode);peak,power=scan_periods(b['t'],b['y'],g,PERIODS,mode)
                    record['models'][mode]=dict(at_gaia=clean(fit),free_peak=peak)
                    write_csv(OUT/f'J0830_K2SC_{extension}_{quality}_{mode}_periodogram.csv',[dict(period_days=float(pi),power=float(vi)) for pi,vi in zip(PERIODS,power)])
                write_csv(OUT/f'J0830_K2SC_{extension}_{quality}_bins.csv',[dict(bjd=float(ti),relative_flux=float(yi),count=int(ci)) for ti,yi,ci in zip(b['t'],b['y'],b['count'])])
                records.append(record)
                if quality=='strict':plots.append(dict(extension=extension,t=b['t'],y=b['y']))
    return records,plots


def downlink_models(curves):
    """Post-review, pipeline-informed offsets; never use target centroids here."""
    records=[];fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for row,target in enumerate(TARGETS):
        data=curves[(target,'TESS','SAP_FLUX','strict')];t,y,g=[data[k] for k in ['t','y','group']];p=TARGETS[target]['period']
        arcs=arc_labels(t,g);record=dict(target=target,arc_count=len(set(arcs)),gap_days=1.5,centroids_or_background_used=False,models={})
        record['segments']=[dict(label=a,first_bjd=float(min(t[arcs==a])),last_bjd=float(max(t[arcs==a])),points=int(sum(arcs==a))) for a in sorted(set(arcs))]
        for mode in ['arc_offset','arc_global_linear','arc_sector_linear','arc_linear']:
            fit=periodic_fit(t,y,g,p,mode);peak,power=scan_periods(t,y,g,PERIODS,mode);holdouts=[]
            for sector in sorted(set(g)):
                train=g!=sector;test=~train;f=periodic_fit(t[train],y[train],g[train],p,mode)
                holdouts.append(dict(sector=sector,training_beta=f['beta'],heldout=heldout_score(t[test],y[test],g[test],p,mode,f['beta'])))
            record['models'][mode]=dict(at_gaia=clean(fit),free_peak=peak,leave_one_sector_out=holdouts)
            write_csv(OUT/f'{target}_TESS_SAP_{mode}_periodogram.csv',[dict(period_days=float(pi),power=float(vi)) for pi,vi in zip(PERIODS,power)])
            axes[row,1].plot(PERIODS,power,label=mode.replace('arc_',''))
            if mode=='arc_sector_linear':
                for sector in sorted(set(g)):
                    m=g==sector;x=t[m]-min(t);axes[row,0].plot(*gap_series(x,(y-fit['_nuisance'])[m]*100,.5),'.',ms=3,label=sector)
                    axes[row,0].plot(*gap_series(x,fit['_periodic'][m]*100,.5),color='black',lw=1)
                write_csv(OUT/f'{target}_TESS_SAP_arc_sector_linear_fit.csv',[dict(bjd=float(ti),sector=str(gi),arc=str(ai),observed_relative_flux=float(yi),fitted_nuisance=float(ni),gaia_period_signal=float(si),residual=float(ri)) for ti,gi,ai,yi,ni,si,ri in zip(t,g,arcs,y,fit['_nuisance'],fit['_periodic'],fit['_residual'])])
        axes[row,0].set(title=f'{target}: six offsets + three sector slopes',xlabel='Days within 2021 epoch',ylabel='Offset/slope-corrected flux (%)');axes[row,0].legend(fontsize=8);axes[row,0].grid(alpha=.2)
        axes[row,1].axvline(p,color='black',ls='--',label=f'Gaia {p:.2f} d');axes[row,1].axvline(p/2,color='gray',ls=':');axes[row,1].set(title=f'{target}: free period scans',xlabel='Period (days)',ylabel='Fractional SSE improvement',xlim=(5,100),ylim=(0,1));axes[row,1].legend(fontsize=7);axes[row,1].grid(alpha=.2)
        records.append(record)
    fig.suptitle('Separate QLP reference offsets at downlink gaps; no centroid/background regressors')
    fig.savefig(OUT/'downlink_offset_models.png',dpi=160);plt.close(fig)
    return records


def main():
    _,curves,_=read_inputs(write_outputs=False)
    robust,corrected=robustness(curves);alternate,kplots=k2sc();arcs=downlink_models(curves)
    result=dict(post_initial_fit=True,reason='Check promising instrument-regression period and ambiguous older K2 variability',regression_robustness=robust,post_review_downlink_models=arcs,downlink_model_status='Added after independent review, motivated by QLP spacecraft-orbit reference images; gap >1.5 d splits each sector into two arcs. Not preregistered.',qlp_reference_documentation='https://doi.org/10.3847/2515-5172/abca2e',k2sc_comparison=alternate,k2sc_documentation='https://archive.stsci.edu/hlsps/k2sc/hlsp_k2sc_k2_llc_all_kepler_v2_readme.txt',limitations=['Conditional models share data and fitted nuisance regressors.','Sector holdout freezes stellar sinusoidal coefficients but refits nuisance terms in the held-out sector.','Target flux-weighted centroids are not independent pointing telemetry.','K2SC is an alternate correction of the same 2015 exposures, not another observing epoch; upstream release 7 and flags differ from the native release 31 product.','EVEREST response exceeded the declared 12 MB retention bound; not retained or analyzed.'])
    (OUT/'followup.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n',encoding='utf8')
    fig,axes=plt.subplots(3,1,figsize=(10,9),layout='constrained')
    for ax,item in zip(axes[:2],corrected):
        for group in sorted(set(item['g'])):
            m=item['g']==group;x=item['t'][m]-min(item['t']);ax.plot(*gap_series(x,item['y'][m]*100,.5),'.',ms=3,label=group);ax.plot(*gap_series(x,item['fit'][m]*100,.5),color='black',lw=1)
        ax.set(title=item['target']+': joint background/centroid model at Gaia period',xlabel='Days within 2021 epoch',ylabel='Model-corrected flux (%)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    for item in kplots:axes[2].plot(*gap_series(item['t']-min(item['t']),item['y']*100,.5),'.-',ms=2,lw=.5,label='K2SC '+item['extension']+', restored time variability')
    axes[2].set(title='J0830: alternate K2 correction, 2015',xlabel='Days within 2015 epoch',ylabel='Relative flux (%)');axes[2].legend(fontsize=8);axes[2].grid(alpha=.2)
    fig.suptitle('Follow-up diagnostics: a jointly fitted correction is a model, not independent confirmation')
    fig.savefig(OUT/'instrument_model_and_k2sc.png',dpi=160);plt.close(fig)
    hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),OUT/'hlsp_k2sc_k2_llc_211569759-c05_kepler_v2_lc.fits'] if p.exists()}
    (OUT/'followup_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf8')
    for r in robust:print(r['target'],'HOLDOUT',[(v['sector'],round(v['heldout']['fraction_improvement'],3)) for v in r['leave_one_sector_out']],flush=True)
    for r in alternate:print('K2SC',r['extension'],r['quality'],{m:round(v['free_peak']['best_period_days'],2) for m,v in r['models'].items()},flush=True)
    for r in arcs:print(r['target'],'DOWNLINK',[(m,round(v['free_peak']['best_period_days'],3),round(v['at_gaia']['harmonic_semiamplitudes'][0]*100,3)) for m,v in r['models'].items()],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--fetch',action='store_true');args=parser.parse_args()
    if args.fetch:fetch_alternatives()
    else:main()
