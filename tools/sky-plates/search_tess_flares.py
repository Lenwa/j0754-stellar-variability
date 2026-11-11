"""Offline, gap-aware exploratory flare search; thresholds select candidates only.

Run from repository root. --injections adds native-time recovery trials.
See the output protocol for predeclared scope and vetting limitations.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analyze_tess_variability as tv

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT/'results/sky-plates-j0754-flare-search-20261110'


def threshold_runs(t, residual, noise):
    t, residual = np.asarray(t), np.asarray(residual)
    if len(t) < 3:
        return []
    cadence = float(np.median(np.diff(t)))
    high = residual > 3*noise
    result = []
    start = None
    for i in range(len(t)+1):
        connected = i < len(t) and high[i]
        if connected and start is not None and t[i]-t[i-1] > 2.5*cadence:
            connected = False
        if start is not None and not connected:
            if i-start >= 3:
                peak = start+int(np.argmax(residual[start:i]))
                if residual[peak] >= 5*noise:
                    result.append(dict(start_index=start, stop_index=i,
                        peak_index=peak, peak_btjd=float(t[peak]),
                        peak_residual=float(residual[peak]),
                        peak_noise_units=float(residual[peak]/noise), points=i-start))
            start = None
        if i < len(t) and high[i] and start is None:
            start = i
    return result


def screen(t, y, width=2.):
    baseline = tv.local_median(t, y, width)
    residual = y-baseline
    dt = np.diff(t)
    short = dt <= 2.5*np.median(dt)
    noise = max(tv.robust_scale(np.diff(residual)[short])/np.sqrt(2), 1e-12)
    return dict(baseline=baseline, residual=residual, noise=noise,
                events=threshold_runs(t, residual, noise),
                negative_events=threshold_runs(t, -residual, noise))


def union_days(groups):
    intervals = sorted(x for group in groups for x in group)
    if not intervals:
        return 0.
    total = 0.; left, right = intervals[0]
    for a, b in intervals[1:]:
        if a <= right+1e-9:
            right = max(right, b)
        else:
            total += right-left
            left, right = a, b
    return float(total+right-left)


def intervals(t):
    half = np.median(np.diff(t))/2
    return [(float(v-half), float(v+half)) for v in t]


def flare_shape(t, peak, fwhm_hours):
    tau = fwhm_hours/24/(.05+np.log(2))
    return tv.event_shape(t, peak, tau, 'flare')


def merged_events(events):
    # Keep the strongest detection in each 0.3-day group, retaining method tags.
    result = []
    for event in sorted(events, key=lambda x:x['peak_noise_units'], reverse=True):
        other = next((r for r in result if abs(r['peak_btjd']-event['peak_btjd']) < .3), None)
        if other is None:
            result.append(dict(event, methods=[event['method']]))
        elif event['method'] not in other['methods']:
            other['methods'].append(event['method'])
    return sorted(result, key=lambda x:x['peak_btjd'])


def contextual_diagnostics(p, event, variants):
    peak = event['peak_btjd']; t = p['b']['t']; y = p['b']['y']
    inside = abs(t-peak) <= .08
    context = abs(t-peak) <= .5
    side = (abs(t-peak) >= .2)&context
    amplitudes = {}
    for width in (1.,2.,3.):
        res = y-tv.local_median(t,y,width)
        amplitudes[str(width)] = float(np.max(res[inside])) if inside.any() else None
    aux = {}
    for name, value in p['ancillary'].items():
        scale = tv.robust_scale(value[side])
        aux[name] = dict(event_minus_side_median=float(np.nanmedian(value[inside])-np.nanmedian(value[side])),
            local_side_MAD_units=float((np.nanmedian(value[inside])-np.nanmedian(value[side]))/scale) if scale and scale > 0 else None,
            flux_correlation=tv.correlation(p['_scan']['residual'][context], (value-tv.local_median(t,value,2.))[context]))
    alternate = {}
    for other in variants:
        tt = other['b']['t']; use = abs(tt-peak) <= .08
        if sum(use) >= 2:
            yy = other['b']['y']; rr = yy-tv.local_median(tt,yy,2.)
            alternate[other['key']] = dict(maximum_residual=float(np.max(rr[use])),
                                          median_residual=float(np.median(rr[use])))
    for name, value in p['apertures'].items():
        rr = value-tv.local_median(t,value,2.)
        alternate[name] = dict(maximum_residual=float(np.max(rr[inside])),median_residual=float(np.median(rr[inside])))
    native = p['native_t']; dt = np.diff(native); cadence=np.median(dt)
    boundary = np.r_[native[0], native[-1], native[:-1][dt>2.5*cadence], native[1:][dt>2.5*cadence]]
    raw = p['raw_t']; q = p['raw_quality']; flagged = np.isfinite(raw)&(q!=0)&(abs(raw-peak)<.25)
    return dict(date_utc=Time(peak+2457000,format='jd',scale='tdb').utc.isot,
        baseline_days_peak_residual=amplitudes,ancillary=aux,alternate_reductions=alternate,
        nearest_quality_gap_or_edge_hours=float(min(abs(boundary-peak))*24),
        flagged_native_rows_within_6h=int(sum(flagged)),
        existing_event=bool(abs(peak-2572.078)<.5),
        assessment='Screened excursion; inspect shape and diagnostics before attribution.')


def search_products(products):
    results=[]
    for p in products:
        if not p['primary']:
            continue
        detections=[]; negatives=[]
        for method,t,y in [('native',p['native_t'],p['native_y']),('20min',p['b']['t'],p['b']['y'])]:
            r=screen(t,y)
            if method=='20min':
                p['_scan']=r
            for sign, target in [('events',detections),('negative_events',negatives)]:
                target.extend(dict(e,method=method,start_btjd=float(t[e['start_index']]),
                    stop_btjd=float(t[e['stop_index']-1])) for e in r[sign])
        events=merged_events(detections)
        for e in events:
            e.update(contextual_diagnostics(p,e,[x for x in products if x['sector']==p['sector'] and x is not p]))
        results.append(dict(key=p['key'],sector=p['sector'],provider=p['provider'],metadata=p['metadata'],
            exposure_days=union_days([intervals(p['native_t'])]),
            noise_20min=p['_scan']['noise'],events=events,negative_events=merged_events(negatives)))
    return results


def injection_trials(products):
    records=[]
    for p in products:
        if not p['primary'] or (p['provider']!='QLP' and p['sector']!=1751):
            continue
        t,y=p['native_t'],p['native_y']; base=screen(t,y); cadence=np.median(np.diff(t))
        b=p['b']; bt=b['t']; bbase=p['_scan']; centers=t.min()+(np.arange(24)+.5)/24*np.ptp(t)
        for width in (.5,2.,4.6):
            tau=width/24/(.05+np.log(2))
            for amplitude in (.005,.01,.025):
                recovered=[]; visible=[]; preexisting=[]
                for peak in centers:
                    shape=flare_shape(t,peak,width); added=amplitude*shape
                    # Retained finite exposures only; do not count gaps as observed.
                    fraction=float(np.sum(shape)*cadence/(tau*(.05+1-np.exp(-6))))
                    visible.append(fraction)
                    rn=tv.injected_residual(t,y,base['baseline'],added,width=2.)
                    bn=tv.bin_values(added,b['inverse'],len(bt))
                    rb=tv.injected_residual(bt,b['y'],bbase['baseline'],bn,width=2.)
                    candidates=threshold_runs(t,rn,base['noise'])+threshold_runs(bt,rb,bbase['noise'])
                    tolerance=max(width/24,1/24)
                    recovered.append(any(abs(c['peak_btjd']-peak)<tolerance for c in candidates))
                    prior=base['events']+bbase['events']
                    preexisting.append(any(abs(c['peak_btjd']-peak)<tolerance for c in prior))
                available=np.asarray(visible)>=.8
                fresh=np.asarray(recovered)&~np.asarray(preexisting)
                records.append(dict(key=p['key'],fwhm_hours=width,amplitude_fractional=amplitude,trials=24,
                    new_recoveries=int(sum(fresh)),preexisting_matches=int(sum(preexisting)),
                    at_least_80percent_profile_sampled=int(sum(available)),
                    new_recoveries_well_sampled=int(sum(fresh&available)),
                    well_sampled_recovery_fraction=float(np.mean(fresh[available])) if available.any() else None,
                    all_time_recovery_fraction=float(np.mean(fresh)),
                    sampled_integral_fraction_p10_p50_p90=np.quantile(visible,[.1,.5,.9]).tolist()))
        print('Injection trials finished:',p['key'],flush=True)
    return records


def plots(products,results,out):
    primary=[p for p in products if p['primary']]
    fig,axes=plt.subplots(len(primary),1,figsize=(13,18),layout='constrained')
    for ax,p,r in zip(axes,primary,results):
        t=p['b']['t']; y=p['_scan']['residual']
        ax.plot(t,y*100,'.',ms=2,alpha=.6)
        for e in r['events']:
            ax.axvline(e['peak_btjd'],color='orange',lw=1)
        ax.axhline(3*r['noise_20min']*100,color='gray',ls=':',lw=.7)
        ax.set(title=r['key'],ylabel='Residual (%)')
    axes[-1].set_xlabel('BTJD; orange: candidate screen, not confirmed flare')
    fig.savefig(out/'all_intervals.png',dpi=130);plt.close(fig)
    events=[(p,e) for p,r in zip(primary,results) for e in r['events']]
    if not events:
        return
    fig,axes=plt.subplots(len(events),2,figsize=(13,max(4,3.2*len(events))),squeeze=False,layout='constrained')
    for (p,e),row in zip(events,axes):
        peak=e['peak_btjd']; use=abs(p['native_t']-peak)<.65; nb=screen(p['native_t'],p['native_y'])
        row[0].plot((p['native_t'][use]-peak)*24,nb['residual'][use]*100,'.',ms=2,alpha=.4,label='Native corrected')
        tt=p['b']['t']; u=abs(tt-peak)<.65
        row[0].plot((tt[u]-peak)*24,p['_scan']['residual'][u]*100,'.-',ms=3,label='20min corrected')
        sap=next(x for x in products if x['file']==p['file'] and not x['primary'])
        st=sap['b']['t']; sy=sap['b']['y']; sr=sy-tv.local_median(st,sy,2.); v=abs(st-peak)<.65
        row[0].plot((st[v]-peak)*24,sr[v]*100,'.-',ms=2,alpha=.7,label='20min SAP')
        row[0].set(title=f"{p['key']} | {e['date_utc'][:10]} | {'known event' if e['existing_event'] else 'new screen candidate'}",ylabel='Residual (%)',xlabel='Hours from detected peak');row[0].legend(fontsize=7)
        for name,values in p['ancillary'].items():
            scale=tv.robust_scale(values[u]); med=np.nanmedian(values[u])
            if scale and scale>0:
                row[1].plot((tt[u]-peak)*24,(values[u]-med)/scale,label=name,lw=.8)
        row[1].set(xlabel='Hours from detected peak',ylabel='Local median/MAD units',ylim=(-12,12),title=f"Nearest quality gap/edge: {e['nearest_quality_gap_or_edge_hours']:.2f}h");row[1].legend(fontsize=6)
    fig.savefig(out/'candidate_diagnostics.png',dpi=130);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=DEFAULT_OUT);ap.add_argument('--injections',action='store_true');args=ap.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    products,inventory=tv.load_products(ROOT/'results/sky-plates-j0754-checks-20261108')
    results=search_products(products)
    primary=[p for p in products if p['primary']]
    coverage=union_days([intervals(p['native_t']) for p in primary])
    preferred=[p for p in primary if p['provider']=='QLP' or p['sector']==1751]
    result=dict(scope='Offline search of existing light curves; no new data retrieved.',inventory=inventory,products=results,
        unique_accepted_cadence_days=coverage,preferred_unique_days=union_days([intervals(p['native_t']) for p in preferred]),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest())
    (out/'search_results.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf8')
    with (out/'all_extracted_screen_hits.csv').open('w',newline='',encoding='utf8') as f:
        columns=['product','sign','peak_btjd','peak_residual','peak_noise_units','methods']
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for r in results:
            for key,sign in [('events','positive'),('negative_events','negative')]:
                for e in r[key]:
                    writer.writerow(dict(product=r['key'],sign=sign,peak_btjd=e['peak_btjd'],
                        peak_residual=e['peak_residual'],peak_noise_units=e['peak_noise_units'],methods=';'.join(e['methods'])))
    plots(products,results,out)
    print(json.dumps(dict(unique_exposure_days=coverage,events=[dict(key=r['key'],events=[dict(peak=e['peak_btjd'],amp=e['peak_residual'],methods=e['methods'],edge_hours=e['nearest_quality_gap_or_edge_hours']) for e in r['events']],negative_count=len(r['negative_events'])) for r in results]),indent=2),flush=True)
    if args.injections:
        trials=injection_trials(products)
        (out/'injections.json').write_text(json.dumps(trials,indent=2,allow_nan=False),encoding='utf8')

if __name__=='__main__':
    main()
