"""Declared quality, ensemble and illustrative wavelength tests for J0754."""
from collections import defaultdict
import math

import numpy as np


def number(value):
    try:
        return float(value)
    except (ValueError,TypeError):
        return float('nan')


def quality(r, band, variant='primary'):
    p=f'w{band}';i=band-1;bad=[]
    def flag(field,expected):
        value=str(r.get(field,''))
        return len(value)>i and value[i]==expected
    for field in ['qual_frame','qi_fact','saa_sep']:
        if not number(r.get(field))>0: bad.append(field)
    if not number(r.get('separation'))<=3: bad.append('separation')
    if not number(r.get('nb'))==1: bad.append('blend')
    for field,expected in [('moon_masked','0'),('cc_flags','0'),('ph_qual','A')]:
        if not flag(field,expected): bad.append(field)
    if not math.isfinite(number(r.get(p+'mpro'))): bad.append('mag')
    error=number(r.get(p+'sigmpro'))
    if not math.isfinite(error) or error<=0: bad.append('error')
    if not number(r.get(p+'snr'))>=10: bad.append('snr')
    if not number(r.get(p+'sat'))==0: bad.append('saturation')
    chi=number(r.get(p+'rchi2'))
    if not (0<=chi<{'primary':3,'loose':5,'strict':2}[variant]): bad.append('chi2')
    if variant!='loose' and not str(r.get(p+'frtr','')).startswith('00'): bad.append('transient')
    return bad


def differential(rows, target, min_controls=3):
    """Exact-frame, unweighted, per-source-centered comparison ensemble.

    No target observation contributes to any comparison-star reference or offset.
    Centering uses all accepted epochs and fixes one offset per comparison star.
    """
    stars=defaultdict(list)
    for r in rows:
        if r['sid']!=target: stars[r['sid']].append(r['mag'])
    centers={s:float(np.median(v)) for s,v in stars.items()}
    frames=defaultdict(dict)
    for r in rows:
        if r['sid']!=target: frames[r['key']][r['sid']]=r['mag']-centers[r['sid']]
    result=[]
    for r in rows:
        if r['sid']!=target: continue
        controls=frames[r['key']]
        if len(controls)<min_controls: continue
        offset=float(np.median(list(controls.values())))
        result.append(dict(r,corrected=r['mag']-offset,offset=offset,n_controls=len(controls)))
    return sorted(result,key=lambda r:r['mjd'])


def visits(rows, field, minimum=5):
    groups=[]
    for r in sorted(rows,key=lambda r:r['mjd']):
        if not groups or r['mjd']-groups[-1][-1]['mjd']>30: groups.append([])
        groups[-1].append(r)
    return [dict(mjd=float(np.median([r['mjd'] for r in g])),
                 first=min(r['mjd'] for r in g),last=max(r['mjd'] for r in g),n=len(g),
                 mag=float(np.median([r[field] for r in g])),
                 scatter=float(1.4826*np.median(np.abs(np.array([r[field] for r in g])-np.median([r[field] for r in g])))))
            for g in groups if len(g)>=minimum]


def sampling_interval(rows,field,start,end,draws=2000):
    """Sampling-only interval: fixed visits, resample complete UTC dates within them.

    This does not estimate visit-to-visit, calibration, comparison-ensemble,
    time-correlated or spectral-model errors. Sparse dates make it discrete.
    """
    groups=[]
    for r in sorted(rows,key=lambda r:r['mjd']):
        if not groups or r['mjd']-groups[-1][-1]['mjd']>30:groups.append([])
        groups[-1].append(r)
    periods=[]
    counts=[]
    for bounds in [start,end]:
        selected=[]
        for g in groups:
            if len(g)<5 or not bounds[0]<=np.median([r['mjd'] for r in g])<bounds[1]:continue
            dates=defaultdict(list)
            for r in g:dates[math.floor(r['mjd'])].append(r[field])
            selected.append(list(dates.values()));counts.append(len(dates))
        periods.append(selected)
    if not all(periods):return dict(status='INSUFFICIENT',dates_per_visit=counts)
    rng=np.random.default_rng(754)
    samples=[]
    for _ in range(draws):
        means=[]
        for period in periods:
            medians=[]
            for date_groups in period:
                choice=rng.integers(0,len(date_groups),len(date_groups))
                medians.append(np.median([x for i in choice for x in date_groups[i]]))
            means.append(np.median(medians))
        samples.append(means[1]-means[0])
    return dict(status='SAMPLING_ONLY',draws=draws,dates_per_visit=counts,
                interval_95=[float(x) for x in np.quantile(samples,[.025,.975])],
                qualification='Fixed visits and ensemble correction; UTC-date resampling within visits only. Excludes calibration, comparison and inter-date correlation errors; not a physical-model confidence interval.')


def extinction_ratio(wavelength,cross_section,at,reference=.48):
    order=np.argsort(wavelength)
    x=np.log(np.asarray(wavelength)[order]);y=np.log(np.asarray(cross_section)[order])
    if min(at,reference)<np.exp(x[0]) or max(at,reference)>np.exp(x[-1]):
        raise ValueError('Extinction wavelength outside tabulated domain')
    return float(np.exp(np.interp(np.log(at),x,y)-np.interp(np.log(reference),x,y)))


def cooling_prediction(temperature, delta_g):
    # Planck B_lambda ratios at fixed wavelength; geometric scale cancels.
    def brightness(w,t): return 1/np.expm1(14387.76877/(w*t))
    def delta(w,t): return -2.5*np.log10(brightness(w,t)/brightness(w,temperature))
    low=temperature/2;high=temperature
    for _ in range(90):
        mid=(low+high)/2
        if delta(.48,mid)>delta_g: low=mid
        else: high=mid
    faint=(low+high)/2
    return dict(temperature_bright=temperature,temperature_faint=faint,
                **{b:float(delta(w,faint)) for b,w in [('g',.48),('r',.62),('W1',3.4),('W2',4.6)]})
