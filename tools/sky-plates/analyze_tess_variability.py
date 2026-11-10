"""Bounded diagnostics of saved J0754 extracted TESS light curves.

Run from repository root: python -B
tools/sky-plates/analyze_tess_variability.py
--out results/sky-plates-j0754-physical-tests-20261108/tess
No downloads, clipping, DC stitching, bolometric flare energies or orbit verdict.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import orbit_photometry as orbit

PERIOD = orbit.PERIOD
BIN = 1/72
DURATIONS = np.array([1, 3, 6])/24


def quality_mask(t, f, q, error=None):
    mask = (np.asarray(q) == 0)&np.isfinite(t)&np.isfinite(f)&(np.asarray(f)>0)
    if error is not None:
        error = np.asarray(error)
        mask &= np.isfinite(error)&(error>0)
    return mask


def robust_scale(y):
    y = np.asarray(y); y = y[np.isfinite(y)]
    return float(1.4826*np.median(abs(y-np.median(y)))) if len(y) else None


def bin_values(y, inverse, n):
    y = np.asarray(y); good = np.isfinite(y)
    total = np.bincount(inverse[good], weights=y[good], minlength=n)
    count = np.bincount(inverse[good], minlength=n)
    return np.divide(total, count, out=np.full(n, np.nan), where=count>0)


def bin_series(t, y, width=BIN):
    key = np.floor((np.asarray(t)-min(t))/width+1e-10).astype(int)
    _, inverse = np.unique(key, return_inverse=True)
    n = max(inverse)+1
    return dict(t=bin_values(t,inverse,n), y=bin_values(y,inverse,n),
                count=np.bincount(inverse), inverse=inverse)


def local_median(t, y, width):
    left = np.searchsorted(t, t-width/2, side='left')
    right = np.searchsorted(t, t+width/2, side='right')
    return np.array([np.median(y[a:b]) for a,b in zip(left,right)])


def equivalent_duration(t, y, max_gap):
    """Signed observed-sampling area in seconds; never bridge missing intervals."""
    dt = np.diff(t); good = (dt>0)&(dt<=max_gap)
    return float(np.sum((np.asarray(y)[1:]+y[:-1])[good]*.5*dt[good])*86400)


def periodic_design(t):
    return orbit.signal(np.asarray(t)+orbit.T0, period=PERIOD, harmonics=2)


def template_fit(t, y, linear=True, groups=None):
    """Joint existing orbit basis with sector nuisance; flux amplitudes, no FAP."""
    t,y = np.asarray(t),np.asarray(y)
    groups = np.zeros(len(t)) if groups is None else np.asarray(groups)
    n = orbit.nuisance(t+orbit.T0, groups, 'group_linear' if linear else 'offset_only')
    s = periodic_design(t)
    # Residualize both data and templates, never detrend only the data first.
    q = np.linalg.svd(n,full_matrices=False)[0][:,:np.linalg.matrix_rank(n)]
    sr = s-q@(q.T@s); yr = y-q@(q.T@y)
    beta = np.linalg.lstsq(sr,yr,rcond=None)[0]
    residual = yr-sr@beta
    retained=[]
    for k in range(2):
        energies=[]
        for phase in np.linspace(0,2*np.pi,24,endpoint=False):
            v=np.array([np.cos(phase),np.sin(phase)])
            energies.append(float(np.sum((sr[:,2*k:2*k+2]@v)**2)/np.sum((s[:,2*k:2*k+2]@v)**2)))
        retained.append([min(energies),max(energies)])
    return dict(beta=beta.tolist(), periods_days=[PERIOD,PERIOD/2],
        fractional_semi_amplitudes=[float(np.linalg.norm(beta[:2])),float(np.linalg.norm(beta[2:]))],
        residual_rms=float(np.sqrt(np.mean(residual**2))),
        null_rms=float(np.sqrt(np.mean(yr**2))), template_condition=float(np.linalg.cond(sr)),
        phase_energy_retained_by_harmonic=retained, linear_sector_nuisance=bool(linear),
        note='Descriptive fixed-template fit; correlated noise and detrending uncertainty not calibrated.')


def box_scores(t, y, noise, durations=DURATIONS):
    """Forward windows on actual occupied bins; reject large internal gaps."""
    t,y = np.asarray(t),np.asarray(y); n=len(t)
    c=np.r_[0,np.cumsum(y)]
    cadence=float(np.median(np.diff(t))) if n>1 else BIN
    gap_limit=max(2.5*cadence,.03)
    gaps=np.r_[0,np.cumsum(np.diff(t)>gap_limit)]
    scores=[]; stops=[]
    for duration in durations:
        end=np.searchsorted(t,t+duration,side='left')
        number=end-np.arange(n); last=np.maximum(0,end-1)
        valid=(number>=3)&((t[last]-t)>=.6*duration)&((gaps[last]-gaps)==0)
        score=np.full(n,np.nan)
        score[valid]=(c[end[valid]]-c[np.flatnonzero(valid)])/(max(float(noise),1e-12)*np.sqrt(number[valid]))
        scores.append(score);stops.append(end)
    return dict(scores=np.array(scores),stops=np.array(stops),durations=np.asarray(durations),gap_limit=gap_limit)


def event_shape(t, onset, duration, kind):
    x=np.asarray(t)-onset
    if kind=='box':
        # Compare absolute endpoints: subtracting onset can round the right
        # endpoint downward and accidentally include a half-open boundary.
        return ((np.asarray(t)>=onset)&(np.asarray(t)<onset+duration)).astype(float)
    result=np.zeros(len(x)); rise=(x>=-.1*duration)&(x<0);decay=(x>=0)&(x<=6*duration)
    result[rise]=1+x[rise]/(.1*duration);result[decay]=np.exp(-x[decay]/duration)
    return result


def bootstrap_screen(t, residual, noise, draws=100, seed=20261108):
    rng=np.random.default_rng(seed); thresholds=[];details=[];n=len(t)
    cadence=np.median(np.diff(t))
    for days in (.25,1.):
        length=max(1,int(round(days/cadence)));maximum=[]
        for _ in range(draws):
            starts=rng.integers(0,n,size=int(np.ceil(n/length)))
            index=np.concatenate([(np.arange(length)+start)%n for start in starts])[:n]
            score=box_scores(t,residual[index],noise)['scores']
            finite=abs(score[np.isfinite(score)]);maximum.append(float(max(finite)) if len(finite) else 0.)
        threshold=float(np.quantile(maximum,.95));thresholds.append(threshold)
        details.append(dict(block_days=days,draws=draws,maximum_score_p50_p95=np.quantile(maximum,[.5,.95]).tolist()))
    return dict(threshold=max(thresholds),blocks=details,
        caveat='Observed-residual block resampling, not a signal-free stationary null or calibrated false alarm probability.')


def correlation(a,b):
    good=np.isfinite(a)&np.isfinite(b)
    if sum(good)<4 or np.std(np.asarray(a)[good])==0 or np.std(np.asarray(b)[good])==0:return None
    return float(np.corrcoef(np.asarray(a)[good],np.asarray(b)[good])[0,1])


def candidates(t,residual,noise,threshold,ancillary):
    scan=box_scores(t,residual,noise);result=[]
    for sign,name in [(1,'positive'),(-1,'negative')]:
        values=sign*scan['scores']; occupied=[]
        order=np.argsort(np.nan_to_num(values,nan=-np.inf).ravel())[::-1]
        for flat in order:
            di,i=np.unravel_index(flat,values.shape)
            if not np.isfinite(values[di,i]):continue
            stop=int(scan['stops'][di,i]);duration=float(scan['durations'][di]);center=(t[i]+t[stop-1])/2
            if any(abs(center-x)<.25 for x in occupied):continue
            occupied.append(center);a=slice(i,stop)
            aux={}
            for k,v in ancillary.items():
                scale=robust_scale(v)
                aux[k]=dict(window_minus_sector_median=float(np.nanmedian(v[a])-np.nanmedian(v)),
                    empirical_MAD_units=float((np.nanmedian(v[a])-np.nanmedian(v))/scale) if scale and scale>0 else None)
            result.append(dict(sign=name,start_btjd=float(t[i]),stop_btjd=float(t[stop-1]),
                duration_hours=duration*24,occupied_bins=int(stop-i),score_white_scaled=float(values[di,i]),
                screening_threshold=threshold,above_screen=bool(values[di,i]>threshold),
                mean_fractional_excursion=float(np.mean(residual[a])),
                signed_observed_equivalent_duration_seconds=equivalent_duration(t[a],residual[a],scan['gap_limit']),
                ancillary_diagnostics=aux,assessment='Unvetted flux excursion, not an identified flare or transit.'))
            if len(occupied)>=3:break
    return result


def injected_residual(t,y,baseline,added,width=1.):
    affected=np.flatnonzero(added!=0)
    if not len(affected):return y-baseline
    changed=(t>=t[affected[0]]-width/2)&(t<=t[affected[-1]]+width/2)
    total=y+added;new=baseline.copy()
    left=np.searchsorted(t,t[changed]-width/2);right=np.searchsorted(t,t[changed]+width/2,side='right')
    new[changed]=[np.median(total[a:b]) for a,b in zip(left,right)]
    return total-new


def injection_grid(native_t,binned,baseline,noise,threshold):
    t,y=binned['t'],binned['y'];inverse=binned['inverse'];n=len(t)
    centers=min(native_t)+(np.arange(24)+.5)/24*np.ptp(native_t)
    cadence=float(np.median(np.diff(native_t)))
    before=box_scores(t,y-baseline,noise)['scores'];records=[]
    for kind,ds in [('box',np.array([1,3,6])/24),('flare',np.array([.5,1,2])/24)]:
        for duration in ds:
            for amp in (.001,.003,.01,.03):
                visible=[];cross=[];pre=[];paired=[];processed=[]
                for onset in centers:
                    shape=event_shape(native_t,onset,duration,kind)
                    integral=duration if kind=='box' else duration*(.05+1-np.exp(-6))
                    visible.append(float(np.sum(shape)*cadence/integral))
                    added=bin_values(shape,inverse,n)*amp*(-1 if kind=='box' else 1)
                    residual=injected_residual(t,y,baseline,added)
                    difference=residual-(y-baseline)
                    processed.append(float(difference@added/(added@added)) if added@added>0 else 0.)
                    scan=box_scores(t,residual,noise)['scores']*(-1 if kind=='box' else 1)
                    # Local search around the true injected event, all declared box widths.
                    local=(t>=onset-max(DURATIONS))&(t<=onset+(duration if kind=='box' else 6*duration))
                    after=scan[:,local];initial=before[:,local]*(-1 if kind=='box' else 1)
                    af=float(np.nanmax(after)) if np.any(np.isfinite(after)) else -np.inf
                    be=float(np.nanmax(initial)) if np.any(np.isfinite(initial)) else -np.inf
                    cross.append(bool(af>threshold));pre.append(bool(be>threshold))
                    paired.append(float(np.max(abs(added))))
                records.append(dict(kind=kind,duration_or_decay_hours=float(duration*24),amplitude_fractional=amp,
                    injected_times=24,threshold_crossings=sum(cross),baseline_crossings=sum(pre),
                    crossing_fraction=float(np.mean(cross)),visible_integral_fraction_p10_p50_p90=np.quantile(visible,[.1,.5,.9]).tolist(),
                    new_threshold_crossings=sum(c and not p for c,p in zip(cross,pre)),
                    new_crossing_fraction=float(np.mean([c and not p for c,p in zip(cross,pre)])),
                    phases_with_under_half_integral_visible=int(sum(v<.5 for v in visible)),
                    binned_peak_response_fraction_p10_p50_p90=np.quantile(np.array(paired)/amp,[.1,.5,.9]).tolist(),
                    paired_processed_projection_gain_p10_p50_p90=np.quantile(processed,[.1,.5,.9]).tolist()))
    return records


def load_products(directory):
    products=[];inventory=[]
    for path in sorted(directory.glob('*18878735*.fits')):
        with fits.open(path) as ff:
            h,d=ff[0].header,ff[1].data;assert int(h['TICID'])==18878735
            provider='QLP' if 'qlp' in path.name else 'TESS-SPOC';sector=int(h['SECTOR'])
            corrected=next(x for x in ('PDCSAP_FLUX','KSPSAP_FLUX','DET_FLUX') if x in d.names)
            inv=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sector=sector,provider=provider,
                columns=[dict(name=c.name,unit=c.unit) for c in ff[1].columns],rows=len(d),quality_counts={str(k):int(v) for k,v in Counter(map(int,d['QUALITY'])).items()},
                data_kind='Extracted aperture/processed lightcurve; no raw telescope pixel time series analyzed.',
                pipeline_version=h.get('PROCVER'),bjdref=ff[1].header.get('BJDREFI',2457000))
            inventory.append(inv)
            for column in ('SAP_FLUX',corrected):
                t=np.asarray(d['TIME'],float);flux=np.asarray(d[column],float)
                error=np.asarray(d[column+'_ERR'],float) if column+'_ERR' in d.names else None
                mask=quality_mask(t,flux,d['QUALITY'],error);order=np.argsort(t[mask]);t0=t[mask][order];f=flux[mask][order];reference=float(np.median(f));norm=f/reference-1
                b=bin_series(t0,norm);aux={}
                apertures={}
                for a in (corrected+'_SML',corrected+'_LAG'):
                    if a in d.names:
                        af=np.asarray(d[a],float)[mask][order]
                        apertures[a]=bin_values(af/np.nanmedian(af)-1,b['inverse'],len(b['t']))
                for a in ('SAP_BKG','SAP_X','SAP_Y','MOM_CENTR1','MOM_CENTR2','POS_CORR1','POS_CORR2'):
                    if a in d.names:aux[a]=bin_values(np.asarray(d[a],float)[mask][order],b['inverse'],len(b['t']))
                products.append(dict(key=f'{provider}_S{sector}_{column}',provider=provider,sector=sector,column=column,
                    file=path.name,primary=column==corrected,native_t=t0,native_y=norm,b=b,ancillary=aux,
                    apertures=apertures,
                    raw_t=t,raw_quality=np.asarray(d['QUALITY']).copy(),
                    metadata=dict(raw_rows=len(d),accepted_rows=int(sum(mask)),rejected_rows=int(sum(~mask)),
                        quality_zero_count=int(sum(d['QUALITY']==0)),native_error_column=column+'_ERR' if error is not None else None,
                        median_relative_native_error=float(np.median(error[mask]/reference)) if error is not None else None,
                        native_normalization_median=reference,first_btjd=float(min(t0)),last_btjd=float(max(t0)),
                        median_cadence_seconds=float(np.median(np.diff(t0))*86400),
                        all_finite_positive_flux_robust_scatter=float(robust_scale(flux[np.isfinite(flux)&(flux>0)]/reference-1)),
                        quality_zero_native_robust_scatter=robust_scale(norm),
                        no_amplitude_clipping=True)))
    return products,inventory


def analyze_product(p,out):
    t,y=p['b']['t'],p['b']['y'];base=local_median(t,y,1.);res=y-base
    noise=max(robust_scale(np.diff(res))/np.sqrt(2),1e-12)
    screen=bootstrap_screen(t,res,noise,seed=int(hashlib.sha256(p['key'].encode()).hexdigest()[:8],16))
    cand=candidates(t,res,noise,screen['threshold'],p['ancillary'])
    for c in cand:
        inside=(p['raw_t']>=c['start_btjd'])&(p['raw_t']<=c['stop_btjd'])
        c['native_nonzero_quality_rows_inside_window']=int(sum(inside&(p['raw_quality']!=0)))
        a=(t>=c['start_btjd'])&(t<=c['stop_btjd'])
        context=(t>=c['start_btjd']-.5)&(t<=c['stop_btjd']+.5)
        c['alternative_aperture_mean_residuals']={k:float(np.mean((v-local_median(t,v,1.))[a])) for k,v in p['apertures'].items()}
        c['local_ancillary_highpass_correlations']={k:correlation(res[context],(v-local_median(t,v,1.))[context]) for k,v in p['ancillary'].items()}
    result=dict(key=p['key'],provider=p['provider'],sector=p['sector'],column=p['column'],primary=p['primary'],file=p['file'],
        metadata=p['metadata'],bins_20min=len(t),binned_scatter=robust_scale(y),
        binned_percentile95_minus5=float(np.diff(np.quantile(y,[.05,.95]))[0]),
        highpass_one_day_scatter=robust_scale(res),highpass_two_day_scatter=robust_scale(y-local_median(t,y,2.)),
        adjacent_difference_noise=noise,lag1_residual_correlation=correlation(res[:-1],res[1:]),
        ancillary_correlations={k:dict(raw=correlation(y,v),highpass=correlation(res,v)) for k,v in p['ancillary'].items()},
        fixed_templates={name:template_fit(t,y,linear) for name,linear in [('constant',False),('linear',True)]},
        screen=screen,candidates=cand)
    fitted=result['fixed_templates']['linear'];s=periodic_design(t);errors=[]
    for harmonic in range(2):
        for phase in np.linspace(0,2*np.pi,12,endpoint=False):
            beta=np.zeros(4);beta[2*harmonic:2*harmonic+2]=.003*np.array([np.cos(phase),np.sin(phase)])
            recovered=np.array(template_fit(t,y+s@beta,True)['beta'])-fitted['beta']
            errors.append(float(np.linalg.norm(recovered-beta)))
    result['paired_fixed_template_injections']=dict(amplitude_fractional=.003,phases_per_harmonic=12,
        maximum_coefficient_response_error=max(errors),
        caveat='Actual retained20min times; tests joint-fit numerics after archive pipeline, not upstream signal transfer or astrophysical coherence.')
    if p['primary']:result['injections']=injection_grid(p['native_t'],p['b'],base,noise,screen['threshold'])
    with (out/(p['key']+'_20min.csv')).open('w',newline='') as f:
        names=['btjd','relative_flux','one_day_baseline','residual','native_samples']+list(p['ancillary']);w=csv.writer(f);w.writerow(names)
        for i in range(len(t)):w.writerow([t[i],y[i],base[i],res[i],p['b']['count'][i]]+[v[i] for v in p['ancillary'].values()])
    p.update(baseline=base,residual=res,result=result)
    return result


def reduction_matches(products):
    matches=[]
    for sector in (71,72):
        a=next(p for p in products if p['sector']==sector and p['provider']=='QLP' and p['primary'])
        b=next(p for p in products if p['sector']==sector and p['provider']=='TESS-SPOC' and p['primary'])
        ta,tb=a['b']['t'],b['b']['t'];idx=np.searchsorted(tb,ta);idx=np.clip(idx,0,len(tb)-1)
        before=np.maximum(0,idx-1);idx=np.where(abs(tb[before]-ta)<abs(tb[idx]-ta),before,idx)
        good=abs(tb[idx]-ta)<BIN*.55
        diff=a['residual'][good]-b['residual'][idx[good]]
        matches.append(dict(sector=sector,matched_bins=int(sum(good)),
            corrected_highpass_correlation=correlation(a['residual'][good],b['residual'][idx[good]]),
            corrected_highpass_difference_robust_scatter=robust_scale(diff),
            caveat='Alternative reductions of same photons; normalization/detrending differences are systematic diagnostics, not independent replication.'))
    return matches


def candidate_concordance(products):
    """Window-level comparisons, preserving alternate reductions as diagnostics."""
    for p in products:
        for c in p['result']['candidates']:
            checks=[]
            for other in products:
                if other is p or other['sector']!=p['sector']:continue
                t=other['b']['t'];inside=(t>=c['start_btjd']-BIN/2)&(t<=c['stop_btjd']+BIN/2)
                if sum(inside)<2:continue
                value=float(np.mean(other['residual'][inside]));same=value*c['mean_fractional_excursion']>0
                checks.append(dict(product=other['key'],bins=int(sum(inside)),mean_residual=value,same_sign=bool(same)))
            c['alternate_reduction_or_SAP_windows']=checks


def plot_with_gaps(ax,t,y,**kwargs):
    t=np.asarray(t);y=np.asarray(y);gap=np.flatnonzero(np.diff(t)>3*BIN)+1
    ax.plot(np.insert(t,gap,np.nan),np.insert(y,gap,np.nan),**kwargs)


def pixel_difference(t,cube,quality,event,target_xy,guard=.15,side=.3):
    """Guarded sideband difference image; descriptive location and motion model."""
    t,cube=np.asarray(t),np.asarray(cube);good=(np.asarray(quality)==0)&np.isfinite(t)&np.all(np.isfinite(cube),axis=(1,2))
    masks=[good&(t>=event[0]-guard-side)&(t<event[0]-guard),
           good&(t>=event[0])&(t<=event[1]),
           good&(t>event[1]+guard)&(t<=event[1]+guard+side)]
    if min(map(np.sum,masks))<3:raise ValueError('Insufficient event/sideband quality frames')
    times=[float(np.mean(t[m])) for m in masks];images=[np.mean(cube[m],axis=0) for m in masks]
    alpha=(times[1]-times[0])/(times[2]-times[0]);baseline=(1-alpha)*images[0]+alpha*images[2];diff=images[1]-baseline
    yy,xx=np.indices(baseline.shape);edge=(xx==0)|(yy==0)|(xx==xx.max())|(yy==yy.max())
    plane=np.column_stack([np.ones(diff.size),xx.ravel(),yy.ravel()]);bg=plane@np.linalg.lstsq(plane[edge.ravel()],diff[edge],rcond=None)[0]
    corrected=diff-bg.reshape(diff.shape);positive=np.maximum(corrected,0)
    centroid=[float(np.sum(xx*positive)/sum(positive.ravel())),float(np.sum(yy*positive)/sum(positive.ravel()))] if sum(positive.ravel())>1e-10 else None
    reference=baseline-float(np.median(baseline[edge]));gy,gx=np.gradient(reference)
    design=np.column_stack([reference.ravel(),gx.ravel(),gy.ravel(),plane])
    coeff=np.linalg.lstsq(design,diff.ravel(),rcond=None)[0];model=(design@coeff).reshape(diff.shape)
    # Compare a stationary-source-only model against the derivative sensitivity.
    fixed=np.column_stack([reference.ravel(),plane]);fc=np.linalg.lstsq(fixed,diff.ravel(),rcond=None)[0]
    core=(xx-target_xy[0])**2+(yy-target_xy[1])**2<=1.5**2
    return dict(quality_counts={name:int(sum(m)) for name,m in zip(['before','event','after'],masks)},
        sideband_mean_btjd=times,interpolation_after_weight=float(alpha),target_wcs_pixel_xy=list(map(float,target_xy)),
        difference_positive_centroid_xy=centroid,
        difference_centroid_offset_pixels=float(np.linalg.norm(np.asarray(centroid)-target_xy)) if centroid else None,
        core_event_fractional_change=float(np.sum(diff[core])/np.sum(reference[core])),
        core_event_fractional_change_definition='Simple aperture contrast: raw event-minus-sideband difference numerator, divided by background-subtracted reference core flux. Not an event-plane-subtracted precision amplitude.',
        stationary_only_fractional_scale=float(fc[0]),
        stationary_plus_motion_model=dict(fractional_stationary_scale=float(coeff[0]),
            derivative_x_coefficient_pixels=float(coeff[1]),derivative_y_coefficient_pixels=float(coeff[2]),
            residual_l2_over_difference_l2=float(np.linalg.norm(diff-model)/max(np.linalg.norm(diff),1e-30))),
        localization_note='Positive-excess moment after outer-ring plane subtraction; no formal position/error distribution or PRF-calibrated astrometry. Derivative model is an instrumental-motion sensitivity, not an orbital motion fit.',
        _maps=dict(before=images[0],event=images[1],after=images[2],baseline=baseline,difference=diff,plane_subtracted_difference=corrected,model_residual=diff-model),
        _masks=masks,_good=good)


def analyze_pixels(path,out):
    event=[2572.0260444438954,2572.1371604654987]
    with fits.open(path,memmap=False) as f:
        d=f[1].data;t=np.asarray(d['TIME'],float);cube=np.asarray(d['FLUX'],float);q=np.asarray(d['QUALITY']);wcs=WCS(f[2].header)
        xy=np.asarray(wcs.all_world2pix([[118.69117410685,16.694621231895]],0)[0]);r=pixel_difference(t,cube,q,event,xy)
        r.update(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sector=f[0].header['SECTOR'],
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            pixel_protocol_sha256=hashlib.sha256((out/'pixel_protocol_amendment.json').read_bytes()).hexdigest() if (out/'pixel_protocol_amendment.json').exists() else None,
            camera=f[0].header['CAMERA'],ccd=f[0].header['CCD'],native_rows=len(t),flux_unit=f[1].columns['FLUX'].unit,
            kind='Native calibrated FFI pixel cutout, not raw unprocessed detector counts.',event_window_btjd=event,
            native_quality_counts={str(k):int(v) for k,v in Counter(map(int,q)).items()},
            available_position_correction_rows={c:int(sum(np.isfinite(d[c]))) for c in ['POS_CORR1','POS_CORR2']})
        yy,xx=np.indices((7,7));core=(xx-xy[0])**2+(yy-xy[1])**2<=1.5**2;edge=(xx==0)|(yy==0)|(xx==6)|(yy==6)
        region=r['_good']&(t>=event[0]-.7)&(t<=event[1]+.7);times=t[region];sub=cube[region]
        border=np.median(sub[:,edge],axis=1);core_flux=np.sum(sub[:,core],axis=1)-border*sum(core.ravel())
        before,eventmask,after=r['_masks'];baseline_flux=np.sum(r['_maps']['baseline'][core])-np.median(r['_maps']['baseline'][edge])*sum(core.ravel())
        with (out/'sector46_pixel_event_timeseries.csv').open('w',newline='') as fp:
            writer=csv.writer(fp);writer.writerow(['btjd','quality','core_background_subtracted_flux_e_per_s','border_median_e_per_s','relative_core_flux','difference_centroid_x','difference_centroid_y','difference_positive_sum_e_per_s','ffi_file'])
            centers=[]
            for i,frame in zip(np.flatnonzero(region),sub):
                a=(t[i]-r['sideband_mean_btjd'][0])/(r['sideband_mean_btjd'][2]-r['sideband_mean_btjd'][0]);base=(1-a)*r['_maps']['before']+a*r['_maps']['after'];diff=frame-base
                plane=np.column_stack([np.ones(49),xx.ravel(),yy.ravel()]);co=plane@np.linalg.lstsq(plane[edge.ravel()],diff[edge],rcond=None)[0];pos=np.maximum(diff-co.reshape(7,7),0);total=float(np.sum(pos));center=[float(np.sum(xx*pos)/total),float(np.sum(yy*pos)/total)] if total else [float('nan')]*2
                k=np.flatnonzero(np.flatnonzero(region)==i)[0];centers.append(center)
                writer.writerow([t[i],int(q[i]),core_flux[k],border[k],core_flux[k]/baseline_flux-1,*center,total,str(d['FFI_FILE'][i])])
        r['local_core_vs_border_correlation']=correlation(core_flux,border)
        # Compare event pixel distribution with target reference distribution.
        ref=r['_maps']['baseline']-np.median(r['_maps']['baseline'][edge]);positive=np.maximum(r['_maps']['plane_subtracted_difference'],0)
        r['reference_positive_centroid_xy']=[float(np.sum(xx*np.maximum(ref,0))/np.sum(np.maximum(ref,0))),float(np.sum(yy*np.maximum(ref,0))/np.sum(np.maximum(ref,0)))]
        r['peak_difference_pixel_yx']=list(map(int,np.unravel_index(np.argmax(positive),positive.shape)))
        r['event_date_utc']=Time(orbit.T0+np.mean(event),format='jd',scale='tdb').utc.isot
        images=r['_maps'];hdus=[fits.PrimaryHDU()]+[fits.ImageHDU(v,header=f[2].header.copy(),name=k.upper()[:16]) for k,v in images.items()];fits.HDUList(hdus).writeto(out/'sector46_event_difference_maps.fits',overwrite=True)
        fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
        for ax,key in zip(axes.ravel()[:4],['baseline','difference','plane_subtracted_difference','model_residual']):
            im=ax.imshow(images[key],origin='lower',cmap='viridis' if key=='baseline' else 'RdBu_r');ax.plot(*xy,'wx',ms=10,mew=2)
            if r['difference_positive_centroid_xy']:ax.plot(*r['difference_positive_centroid_xy'],'k+',ms=10,mew=2)
            ax.set(title=key.replace('_',' '),xlabel='Cutout x pixel',ylabel='Cutout y pixel');fig.colorbar(im,ax=ax,label='e-/s')
        ax=axes.ravel()[4];ax.plot((times-np.mean(event))*24,core_flux/baseline_flux-1,'.-');ax.axvspan((event[0]-np.mean(event))*24,(event[1]-np.mean(event))*24,color='orange',alpha=.2);ax.set(xlabel='Hours from selected event center',ylabel='Core aperture relative flux',title='New pixel extraction; quality=0')
        ax=axes.ravel()[5];centers=np.array(centers);im=ax.scatter(centers[:,0],centers[:,1],c=(times-np.mean(event))*24,s=12);ax.plot(*xy,'rx',ms=10);ax.set(xlabel='Positive difference x',ylabel='Positive difference y',title='All local difference moments (noise included)');fig.colorbar(im,ax=ax,label='Hours from event center')
        fig.suptitle('Sector46 calibrated7x7pixels: Gaia WCS white cross; positive-excess centroid black plus')
        fig.savefig(out/'sector46_pixel_localization.png',dpi=150);plt.close(fig)
    clean={k:v for k,v in r.items() if not k.startswith('_')};(out/'sector46_pixel_assessment.json').write_text(json.dumps(clean,indent=2,allow_nan=False),encoding='utf8');return clean


def figures(products,out):
    main=[p for p in products if p['primary']]
    fig,axes=plt.subplots(len(main),1,figsize=(13,20),layout='constrained')
    for ax,p in zip(axes,main):
        t,y=p['b']['t'],p['b']['y'];sap=next(x for x in products if x['file']==p['file'] and not x['primary'])
        ax.plot(sap['b']['t'],sap['b']['y'],'.',ms=1,alpha=.3,label='SAP separately normalized')
        ax.plot(t,y,'.',ms=2,alpha=.5,label=p['column']);plot_with_gaps(ax,t,p['baseline'],lw=1,label='1-day running median')
        ax.set(title=p['key'],ylabel='Relative flux');ax.legend(fontsize=7,ncol=3)
    axes[-1].set_xlabel('BTJD = BJD_TDB - 2457000; gaps remain open')
    fig.savefig(out/'within_sector_variability.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    for ax,p in zip(axes.ravel(),main):
        for kind,style in [('box','-'),('flare','--')]:
            rr=[r for r in p['result']['injections'] if r['kind']==kind]
            for duration in sorted(set(r['duration_or_decay_hours'] for r in rr)):
                a=[r for r in rr if r['duration_or_decay_hours']==duration]
                ax.semilogx([r['amplitude_fractional'] for r in a],[r['new_crossing_fraction'] for r in a],style,marker='.',label=f'{kind} {duration:g}h')
        ax.set(title=p['key'],ylim=(-.05,1.05),xlabel='Injected fractional amplitude',ylabel='Screen crossings / 24 times');ax.legend(fontsize=6)
    fig.suptitle('Actual native-time injections before 20-minute binning; crossings are not universal completeness limits')
    fig.savefig(out/'event_injection_response.png',dpi=150);plt.close(fig)
    top=sorted([(c['score_white_scaled']/p['result']['screen']['threshold'],p,c) for p in main for c in p['result']['candidates']],key=lambda v:v[0],reverse=True)[:6]
    fig,axes=plt.subplots(3,2,figsize=(14,11),layout='constrained')
    for ax,(_,p,c) in zip(axes.ravel(),top):
        t=p['b']['t'];center=(c['start_btjd']+c['stop_btjd'])/2;use=abs(t-center)<.6
        # Gap plotting uses days before converting its x coordinates to hours.
        tt=t[use]-center;yy=p['residual'][use];g=np.flatnonzero(np.diff(tt)>3*BIN)+1
        ax.plot(np.insert(tt,g,np.nan)*24,np.insert(yy,g,np.nan),'.-',ms=3,lw=.5)
        ax.axvspan((c['start_btjd']-center)*24,(c['stop_btjd']-center)*24,color='orange',alpha=.2)
        ax.set(title=f"{p['key']} {c['sign']} screen ratio={c['score_white_scaled']/c['screening_threshold']:.2f}",xlabel='Hours from candidate center',ylabel='Relative flux minus 1-day median')
    fig.savefig(out/'top_flux_excursions.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=Path('results/sky-plates-j0754-checks-20261108'));parser.add_argument('--out',type=Path,required=True);parser.add_argument('--pixels-only',type=Path);args=parser.parse_args();out=args.out
    if not (out/'protocol.json').exists():raise RuntimeError('Save protocol before analysis')
    if args.pixels_only:
        if not (out/'pixel_protocol_amendment.json').exists():raise RuntimeError('Save pixel protocol before inspection')
        print(json.dumps(analyze_pixels(args.pixels_only,out),indent=2));return
    products,inventory=load_products(args.input);results=[]
    for p in products:
        print('Analyzing',p['key'],flush=True);results.append(analyze_product(p,out))
    candidate_concordance(products)
    epoch_templates={}
    for column in ('SAP_FLUX','corrected'):
        for year,sectors in [('2021',(44,45,46)),('2023',(71,72))]:
            pp=[p for p in products if p['provider']=='QLP' and p['sector'] in sectors and ((p['primary']) if column=='corrected' else p['column']==column)]
            t=np.concatenate([p['b']['t'] for p in pp]);y=np.concatenate([p['b']['y'] for p in pp]);groups=np.concatenate([np.full(len(p['b']['t']),p['sector']) for p in pp])
            epoch_templates[column+'_'+year]={name:template_fit(t,y,linear,groups) for name,linear in [('constant',False),('linear',True)]}
    result=dict(inventory=inventory,products=results,reduction_matches=reduction_matches(products),epoch_templates=epoch_templates,
        prior_fixed_period_result='results/sky-plates-j0754-models-20261108/short_period/report.md',
        scope='New within-sector and event-sensitivity diagnostics on existing extracted products; not raw-pixel verification or a new orbital rejection.',
        protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf8');figures(products,out)
    report(result,out)


def report(r,out):
    lines=['J0754 TESS extracted-lightcurve physical tests','',
        'OUTCOME: The archive products show strong reduction dependence and contain an unvetted',
        'flare-like positive excursion in Sector 46. They do not validate a universal 0.002 mag',
        'transit limit or a bolometric flare-energy threshold. The screen below resamples its',
        'own candidates: its threshold often equals the largest observed excursion. No events',
        'above that self-including descriptive screen is NOT evidence of flare absence.','',
        'These eight FITS files are extracted light curves, not raw telescope pixel data.',
        'QUALITY=0 and finite positive flux/error screens were used without brightness clipping.',
        'Each sector/product was normalized separately; no multi-year DC level is inferred.','',
        'Primary corrected products; empirical fractional-flux scatter, not a precision limit:',
        'Product | native accepted | 20min bins | robust scatter | 1day highpass scatter | lag1 correlation']
    for p in r['products']:
        if p['primary']:lines.append(f"{p['key']} | {p['metadata']['accepted_rows']} | {p['bins_20min']} | {p['binned_scatter']:.5g} | {p['highpass_one_day_scatter']:.5g} | {p['lag1_residual_correlation']:.3f}")
    lines+=['','Background and centroid Pearson correlations, two-day baseline sensitivity, native errors,',
        'quality counts and constant/linear fixed-template fits are retained in results.json.',
        'Correlations do not by themselves prove contamination or astrophysical variability.','',
        'Flux excursions: top three of each sign per product are retained regardless of screening.',
        'Their white-scaled scores are not Gaussian detection significances. The screening threshold',
        'is the larger 95th percentile of maximum statistics from 6 h/24 h observed-residual block',
        'resampling. Nonstationarity, candidate signals in the resampling pool and residual systematics',
        'prevent treating this as a calibrated null false alarm probability.','']
    for p in r['products']:
        if p['primary']:
            above=[c for c in p['candidates'] if c['above_screen']]
            lines.append(f"{p['key']}: {len(above)} of 6 retained excursions above descriptive screen; threshold={p['screen']['threshold']:.3g}.")
    lines+=['','The injection grid uses actual quality-retained native timestamps, injects before 20 minute',
        'binning, and repeats the 1 day local median. Box deficits 1/3/6 h and flare-like decay 0.5/1/2 h',
        'at fractional amplitudes 0.001/0.003/0.01/0.03 are tested at 24 times per product, including',
        'gaps. Threshold-crossing fractions, baseline crossings and visible integral fractions are',
        'recorded separately. These are bounded sampling/pipeline sensitivities, not universal',
        'upper limits or a correction for upstream pipeline attenuation.','',
        'The 1/3/6 hour grid concerns brief events. Multi-day giant-binary transits can last',
        'about 2.7â€“5.7 days in the separate conditional models. The 1 day highpass suppresses such',
        'broad events, so this grid supplies NO transit-depth/radius bound for those geometries.','',
        'A fixed 35.9204566 day fundamental and 17.9602283 day harmonic are jointly fitted with',
        'constant or linear sector nuisances using the existing orbit_photometry basis. Short',
        'baselines and small retained signal energy can make amplitudes poorly identified. Fits',
        'and reduction comparisons do not establish a stable waveform, rotation or a new orbit',
        'exclusion. The existing independent ASAS-SN/ZTF period result remains separate.','',
        'Equivalent durations integrate signed observed residuals in seconds and never bridge',
        'large gaps. They are window-dependent extracted-aperture quantities. Bolometric flare',
        'energies require an unknown flare/passband correction, quiescent band luminosity and',
        'dilution model. No 150 Lsun, 0.002 mag limit, 100 ppm precision or universal 3 sigma bound was',
        'assumed. Transit radius limits require validated completeness and geometry; none is claimed.','',
        'Files: results.json; per-product 20 min CSV; within_sector_variability.png;',
        'event_injection_response.png; top_flux_excursions.png. Input/source/protocol hashes are',
        'retained. No new FFI downloads, telescope actions, accounts or publication occurred.']
    lines[-1]='retained. No complete FFI files, telescope actions, accounts or publication occurred.'
    lines+=['','OVERLAPPING CORRECTED REDUCTIONS (SAME PHOTONS)']
    for item in r['reduction_matches']:
        lines.append(f"Sector {item['sector']}: highpass correlation {item['corrected_highpass_correlation']:.3f}; robust difference scatter {item['corrected_highpass_difference_robust_scatter']:.5g}.")
    lines+=['','FIXED-TEMPLATE FRACTIONAL SEMI-AMPLITUDES; DESCRIPTIVE, NOT UPPER LIMITS']
    for name,variants in r['epoch_templates'].items():
        for mode,fit in variants.items():lines.append(f"{name} {mode}: P={fit['fractional_semi_amplitudes'][0]:.5g}, P/2={fit['fractional_semi_amplitudes'][1]:.5g}; condition={fit['template_condition']:.3g}.")
    if (out/'sector46_pixel_assessment.json').exists():
        pixel=json.loads((out/'sector46_pixel_assessment.json').read_text())
        lines+=['','BOUNDED SECTOR 46 PIXEL FOLLOW-UP (SEPARATE PROTOCOL AMENDMENT)',
            'A 3.91 MB official TESSCut 7x7 calibrated FFI cutout was subsequently recovered.',
            f"Quality-zero frames before/event/after: {pixel['quality_counts']}; event date {pixel['event_date_utc']}.",
            f"Positive difference centroid offset from Gaia WCS: {pixel['difference_centroid_offset_pixels']:.3f}pixel; no calibrated positional confidence interval.",
            f"Stationary-source-shape scale: {pixel['stationary_only_fractional_scale']:.5f}; with derivative motion terms: {pixel['stationary_plus_motion_model']['fractional_stationary_scale']:.5f}.",
            'The excursion spans many native frames and is spatially compact near the target.',
            'This improves its status to a target-associated flare-like candidate, but does not',
            'exclude unresolved blends, a near-aligned moving source, PRF/background effects or',
            'calibration artifacts. See pixel maps, per-frame table and sideband sensitivities.',
            'It is not a confirmed stellar flare or a measurement of bolometric flare energy.']
    (out/'report.txt').write_text('\n'.join(lines)+'\n',encoding='utf8')


if __name__=='__main__':main()
