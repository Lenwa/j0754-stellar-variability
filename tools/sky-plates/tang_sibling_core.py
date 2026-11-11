"""Numerical tools for short, separately calibrated stellar epochs."""
from collections import defaultdict
import numpy as np


def absolute_time(relative,header):
    return np.asarray(relative,float)+float(header['BJDREFI'])+float(header.get('BJDREFF',0))+float(header.get('TIMEZERO',0))


def quality_mask(flags,mode='strict'):
    flags=np.asarray(flags,np.int64)
    if mode=='strict':return flags==0
    if mode=='kepler_default':return (flags&1130799)==0
    raise ValueError('Unknown quality selection')


def restore_k2sc_variability(flux,time_trend):
    return np.asarray(flux,float)+np.asarray(time_trend,float)-np.nanmedian(time_trend)


def bin_series(t,y,groups,width=.2,extra=None):
    """Collapse exact duplicates then median-bin within each calibration group."""
    t,y,groups=np.asarray(t,float),np.asarray(y,float),np.asarray(groups)
    extra=np.empty((len(t),0)) if extra is None else np.asarray(extra,float)
    exact=defaultdict(list)
    for i,(ti,yi,g) in enumerate(zip(t,y,groups)):
        if np.isfinite(ti) and np.isfinite(yi):exact[(str(g),float(ti))].append(i)
    bins=defaultdict(list)
    for (g,ti),indices in exact.items():
        values=np.r_[ti,np.median(y[indices]),np.nanmedian(extra[indices],axis=0)]
        bins[(g,int(np.floor(ti/width)))].append(values)
    rows=[(g,np.median(values,axis=0),len(values)) for (g,_),values in sorted(bins.items(),key=lambda item:(item[0][0],item[0][1]))]
    if not rows:raise ValueError('No valid data to bin')
    return dict(t=np.array([v[0] for _,v,_ in rows]),y=np.array([v[1] for _,v,_ in rows]),group=np.array([g for g,_,_ in rows]),extra=np.array([v[2:] for _,v,_ in rows]),count=np.array([n for _,_,n in rows]))


def orthogonal_columns(matrix):
    u,s,_=np.linalg.svd(matrix,full_matrices=False)
    rank=int(np.sum(s>max(matrix.shape)*np.finfo(float).eps*s[0])) if len(s) and s[0]>0 else 0
    return u[:,:rank]


def arc_labels(t,groups,gap_days=1.5):
    """Split each sector at observed downlink gaps; preserve input order.

    This is a time-based proxy for separate QLP spacecraft-orbit references,
    not a claim that each short quality gap defines a new photometric reference.
    """
    t,groups=np.asarray(t,float),np.asarray(groups);labels=np.empty(len(t),object)
    for group in np.unique(groups):
        indices=np.flatnonzero(groups==group);indices=indices[np.argsort(t[indices])]
        arcs=np.r_[0,np.cumsum(np.diff(t[indices])>gap_days)]
        labels[indices]=[f'{group}_arc{int(arc)}' for arc in arcs]
    return labels.astype(str)


def nuisance_matrix(t,groups,mode,extra=None):
    t,groups=np.asarray(t,float),np.asarray(groups)
    columns=[]
    if mode in ['arc_offset','arc_global_linear','arc_sector_linear','arc_linear']:
        arcs=arc_labels(t,groups)
        for arc in np.unique(arcs):
            mask=arcs==arc;columns.append(mask.astype(float))
            if mode=='arc_linear':columns.append(mask*(t-np.mean(t[mask]))/max(1,np.ptp(t[mask])))
        if mode=='arc_global_linear':columns.append((t-np.mean(t))/max(1,np.ptp(t)))
        if mode=='arc_sector_linear':
            for group in np.unique(groups):
                mask=groups==group;columns.append(mask*(t-np.mean(t[mask]))/max(1,np.ptp(t[mask])))
        return np.column_stack(columns)
    for group in np.unique(groups):
        mask=groups==group;x=(t-np.mean(t[mask]))/max(1,np.ptp(t[mask]))
        columns.append(mask.astype(float))
        if mode in ['sector_linear','sector_quadratic','instrument']:
            columns.append(mask*x)
        if mode=='sector_quadratic':columns.append(mask*x*x)
        if mode=='instrument' and extra is not None:
            for col in np.asarray(extra).T:
                values=col[mask]
                if not np.all(np.isfinite(values)):raise ValueError('Instrument regressors must be finite')
                scatter=np.std(values)
                if scatter>0:columns.append(mask*(col-np.mean(values))/scatter)
    if mode=='global_linear':columns.append((t-np.mean(t))/max(1,np.ptp(t)))
    if mode not in ['offset','global_linear','sector_linear','sector_quadratic','instrument']:raise ValueError('Unknown nuisance model')
    return np.column_stack(columns)


def signal_matrix(t,period,harmonics=1):
    angle=2*np.pi*(np.asarray(t,float)-2457000.)/period
    return np.column_stack([fun(k*angle) for k in range(1,harmonics+1) for fun in [np.cos,np.sin]])


def periodic_fit(t,y,groups,period,mode='global_linear',harmonics=1,extra=None):
    t,y=np.asarray(t,float),np.asarray(y,float)
    n=nuisance_matrix(t,groups,mode,extra);q=orthogonal_columns(n)
    s=signal_matrix(t,period,harmonics);sr=s-q@(q.T@s);yr=y-q@(q.T@y)
    beta=np.linalg.lstsq(sr,yr,rcond=None)[0]
    residual=yr-sr@beta;null=float(yr@yr)
    scale=max(float(y@y),1.)
    power=1-float(residual@residual)/null if null>1e-24*scale else 0.
    eigen,vectors=np.linalg.eigh(s.T@s)
    invroot=(vectors/np.sqrt(np.maximum(eigen,1e-30)))@vectors.T
    retention=np.clip(np.linalg.eigvalsh(invroot@(sr.T@sr)@invroot),0,1)
    return dict(period_days=float(period),harmonics=harmonics,nuisance=mode,nuisance_rank=q.shape[1],beta=beta.tolist(),harmonic_semiamplitudes=[float(np.linalg.norm(beta[2*k:2*k+2])) for k in range(harmonics)],power=float(power),null_sse=null,residual_sse=float(residual@residual),phase_energy_retention=retention.tolist(),span_cycles=float(np.ptp(t)/period),signal_condition=float(np.linalg.cond(sr)),_periodic=s@beta,_nuisance=y-residual-s@beta,_residual=residual)


def scan_periods(t,y,groups,periods,mode='global_linear',harmonics=1,extra=None,return_bases=False):
    q=orthogonal_columns(nuisance_matrix(t,groups,mode,extra));yr=y-q@(q.T@y);denom=float(yr@yr)
    basis=[]
    for period in periods:
        a=signal_matrix(t,period,harmonics);a-=q@(q.T@a)
        basis.append(orthogonal_columns(a).T)
    basis=np.array(basis)
    flat=denom<=1e-24*max(float(np.asarray(y)@np.asarray(y)),1.)
    power=np.zeros(len(periods)) if flat else np.sum((basis@yr)**2,axis=1)/denom
    index=int(np.argmax(power))
    result=dict(best_period_days=None if flat else float(periods[index]),power=float(power[index]),span_cycles=None if flat else float(np.ptp(t)/periods[index]),at_search_edge=False if flat else bool(index in [0,len(periods)-1]),flat_null=bool(flat))
    if return_bases:return result,power,basis,q,yr
    return result,power


def heldout_score(t,y,groups,period,mode,beta,extra=None):
    q=orthogonal_columns(nuisance_matrix(t,groups,mode,extra))
    observed=y-q@(q.T@y)
    s=signal_matrix(t,period,len(beta)//2);s-=q@(q.T@s)
    residual=observed-s@np.asarray(beta);null=float(observed@observed)
    return dict(fraction_improvement=float(1-residual@residual/null) if null>1e-24 else 0.,null_sse=null,frozen_sse=float(residual@residual),fit_only_nuisance_on_heldout=True,retention=periodic_fit(t,y,groups,period,mode,len(beta)//2,extra)['phase_energy_retention'])


def count_exceedances(values,observed):
    """Include algebraic ties despite projection/BLAS round-off."""
    return int(np.sum(np.asarray(values)>=observed-1e-12*max(1,abs(observed))))


def block_sign_diagnostic(t,y,groups,periods,mode='global_linear',extra=None,draws=256,block_days=(7,14),seed=20261111):
    result,power,basis,q,yr=scan_periods(t,y,groups,periods,mode,extra=extra,return_bases=True)
    rng=np.random.default_rng(seed);outputs={}
    for days in block_days:
        block=np.floor((t-min(t))/days).astype(int)
        signs=rng.choice([-1.,1.],size=(max(block)+1,draws))
        z=yr[:,None]*signs[block];z-=q@(q.T@z)
        norm=np.sum(z*z,axis=0)
        simulated=np.sum((basis@z)**2,axis=1)/np.maximum(norm,1e-300)[None,:]
        count=count_exceedances(np.max(simulated,axis=0),result['power'])
        entry=dict(draws=draws,occupied_blocks=int(len(set(block))),trial_exceedances=count,raw_trial_fraction=count/draws,fraction_scanned_max_ge_observed=float((1+count)/(draws+1)),monte_carlo_convention='Plus-one correction; numerical ties counted at relative/absolute tolerance 1e-12',status='Conditional block-sign sensitivity, not a calibrated false-alarm probability')
        unique,inverse=np.unique(block,return_inverse=True)
        if len(unique)<=8:
            patterns=((np.arange(2**len(unique))[:,None]>>np.arange(len(unique)))&1)*2-1
            exact=yr[:,None]*patterns[:,inverse].T;exact-=q@(q.T@exact)
            exact_norm=np.sum(exact*exact,axis=0)
            exact_max=np.max(np.sum((basis@exact)**2,axis=1)/np.maximum(exact_norm,1e-300)[None,:],axis=0)
            n=count_exceedances(exact_max,result['power'])
            entry.update(exact_patterns=len(patterns),exact_exceedances=n,exact_conditional_fraction=n/len(patterns))
        outputs[str(days)]=entry
    return outputs


def local_residual_scale(t,residual,half_window=.5,minimum=20):
    t,residual=np.asarray(t,float),np.asarray(residual,float)
    scales=[]
    for ti in t:
        values=residual[(np.abs(t-ti)<=half_window)&np.isfinite(residual)]
        scale=1.4826*np.median(np.abs(values-np.median(values))) if len(values)>=minimum else np.nan
        scales.append(max(scale,1e-12) if np.isfinite(scale) else np.nan)
    return np.asarray(scales)


def excursion_runs(t,above,cadence,minimum=3,cadence_numbers=None,cadence_step=1):
    result=[];run=[]
    for i,(ti,hit) in enumerate(zip(t,above)):
        gap=bool(run and ti-t[run[-1]]>2.5*cadence)
        if run and cadence_numbers is not None:gap|=cadence_numbers[i]-cadence_numbers[run[-1]]!=cadence_step
        if not hit or gap:
            if len(run)>=minimum:result.append(run)
            run=[]
        if hit:run.append(i)
    if len(run)>=minimum:result.append(run)
    return result


def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items() if not k.startswith('_')}
    if isinstance(value,list):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    return value
