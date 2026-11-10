"""Reproduce the predeclared J0754 optical/infrared diagnostic; no network."""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from asassn_replication import parse_row
from dust_models import quality,differential,visits,extinction_ratio,cooling_prediction,sampling_interval

TARGET='667401913089980544'


def stamp(v):
    d=Time(v['mjd'],format='mjd').to_datetime()
    return dict(v,year=d.year,month=d.month,half=1 if d.month<=6 else 2,
                date=d.isoformat(),decimal_year=float(Time(v['mjd'],format='mjd').decimalyear))


def contrast(v,start=2018,end=2022):
    a=[x['mag'] for x in v if x['year']==start]
    b=[x['mag'] for x in v if x['year']==end]
    if not a or not b: return dict(status='INSUFFICIENT',visits_start=len(a),visits_end=len(b))
    half=[]
    for h in [1,2]:
        x=[r['mag'] for r in v if r['year']==start and r['half']==h]
        y=[r['mag'] for r in v if r['year']==end and r['half']==h]
        if x and y: half.append(dict(half=h,delta=float(np.median(y)-np.median(x))))
    return dict(status='MEASURED',delta=float(np.median(b)-np.median(a)),
                visits_start=len(a),visits_end=len(b),matched_halves=half,
                seasonal_support='MATCHED' if half else 'NO_COMMON_HALF_YEAR')


def main(previous,asassn,out):
    gaia={r['source_id']:r for r in csv.DictReader((previous/'gaia_field.csv').open())}
    ids=json.loads((previous/'isotonic_pixel_results.json').read_text())['results'][-1]['selected_controls']+[TARGET]
    native=[];dedup=[];counts={};hashes={}
    for sid in ids:
        path=out/f'neowise_{sid}.csv';data=list(csv.DictReader(path.open()))
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        r0=gaia[sid];c0=SkyCoord(float(r0['ra'])*u.deg,float(r0['dec'])*u.deg)
        grouped=defaultdict(list)
        for r in data:
            r=dict(r,sid=sid,key=(r['scan_id'],r['frame_num']),mjd=float(r['mjd']))
            r['separation']=float(c0.separation(SkyCoord(float(r['ra'])*u.deg,float(r['dec'])*u.deg)).arcsec)
            grouped[r['key']].append(r);native.append(r)
        # Astrometric choice only, before any photometric quality or brightness use.
        selected=[min(g,key=lambda r:r['separation']) for g in grouped.values()]
        dedup.extend(selected)
        counts[sid]=dict(raw=len(data),frames=len(selected),
                         min_separation=min([r['separation'] for r in selected],default=None),
                         max_separation=max([r['separation'] for r in selected],default=None))
    results={};plot_data={}
    for variant in ['primary','strict','loose']:
        bands={};target_raw={};target_corrected={}
        for band in [1,2]:
            accepted=[];rejected=defaultdict(Counter)
            for r in dedup:
                bad=quality(r,band,variant)
                if bad: rejected[r['sid']].update(bad)
                else: accepted.append(dict(r,mag=float(r[f'w{band}mpro'])))
            raw=[r for r in accepted if r['sid']==TARGET]
            corrected=differential(accepted,TARGET)
            vr=[stamp(v) for v in visits(raw,'mag')]
            vc=[stamp(v) for v in visits(corrected,'corrected')]
            sample=sampling_interval(corrected,'corrected',
                (Time('2018-01-01').mjd,Time('2019-01-01').mjd),
                (Time('2022-01-01').mjd,Time('2023-01-01').mjd))
            placebos={}
            for sid in ids:
                if sid==TARGET: continue
                p=differential([r for r in accepted if r['sid']!=TARGET],sid)
                placebos[sid]=contrast([stamp(v) for v in visits(p,'corrected')])
            jackknife=[]
            for sid in ids:
                if sid==TARGET:continue
                c=differential([r for r in accepted if r['sid']!=sid],TARGET)
                z=contrast([stamp(v) for v in visits(c,'corrected')])
                if z['status']=='MEASURED': jackknife.append(dict(omitted=sid,delta=z['delta']))
            bands[f'W{band}']=dict(accepted_by_star=dict(Counter(r['sid'] for r in accepted)),
                rejected_reasons={s:dict(v) for s,v in rejected.items()},raw_visits=vr,corrected_visits=vc,
                raw_contrast=contrast(vr),corrected_contrast=contrast(vc),
                raw_recovery=contrast(vr,2022,2024),corrected_recovery=contrast(vc,2022,2024),
                target_raw_frames=len(raw),target_corrected_frames=len(corrected),
                sampling_only_interval=sample,
                comparison_counts=dict(Counter(r['n_controls'] for r in corrected)),
                control_placebos=placebos,comparison_jackknife=jackknife)
            target_raw[band]=raw;target_corrected[band]=corrected
        # Color is computed only from paired bands on the same native exposure.
        colors={}
        for label,streams,field in [('raw',target_raw,'mag'),('corrected',target_corrected,'corrected')]:
            b={r['key']:r for r in streams[2]}
            paired=[dict(mjd=r['mjd'],mag=r[field]-b[r['key']][field])
                    for r in streams[1] if r['key'] in b]
            v=[stamp(x) for x in visits(paired,'mag')]
            colors[label]=dict(frames=len(paired),visits=v,contrast=contrast(v))
        results[variant]=dict(bands=bands,colors=colors)
        if variant=='primary':plot_data=target_corrected
    models={}
    for name in ['rv31','rv4','rv5.5']:
        p=out/f'draine_{name}_extinction.txt';grid=[]
        for line in p.read_text().splitlines():
            try:
                values=[float(v) for v in line.split()]
                if len(values)==6 and values[0]>0 and values[3]>0: grid.append(values)
            except ValueError: pass
        grid=np.array(grid)
        models[name]={b:.30*extinction_ratio(grid[:,0],grid[:,3],w) for b,w in [('g',.48),('r',.62),('W1',3.4),('W2',4.6)]}
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    models['fixed_area_cooling']=[cooling_prediction(t,.30) for t in [4000,4500,5000,5500]]
    optical=[]
    for row in csv.DictReader(asassn.open()):
        r=parse_row(row)
        if r['accepted'] and r['band']=='g':optical.append(r)
    optical_matches={}
    irvis=results['primary']['bands']['W1']['corrected_visits']
    for window in [30,15]:
        matched=[]
        for v in irvis:
            jd=v['mjd']+2400000.5
            for camera in sorted({r['camera'] for r in optical}):
                g=[r for r in optical if r['camera']==camera and abs(r['jd']-jd)<=window]
                dates=len({r['day'] for r in g})
                if dates>=3:
                    matched.append(dict(v,camera=camera,mag=float(np.median([r['mag'] for r in g])),
                                        infrared_mag=v['mag'],
                                        optical_dates=dates,optical_points=len(g)))
        cameras=sorted({r['camera'] for r in matched})
        optical_matches[str(window)]=dict(matches=matched,
            contrasts={c:contrast([v for v in matched if v['camera']==c]) for c in cameras},
            g_minus_W1_changes={c:contrast([dict(v,mag=v['mag']-v['infrared_mag']) for v in matched if v['camera']==c]) for c in cameras})
    hashes[str(asassn)]=hashlib.sha256(asassn.read_bytes()).hexdigest()
    report=dict(target=TARGET,native_counts=counts,variants=results,models=models,optical_matches=optical_matches,
        limitations=['Visit medians are not independent stellar-state replicates; two seasonal visits/year give no reliable total systematic uncertainty.',
                     'Exact-frame corrections assume the fixed comparison ensemble has a stable median; placebos and leave-one-control-out sensitivities are diagnostics, not proof.',
                     'Monochromatic extinction and fixed-area Planck models are illustrations, not stellar atmosphere or circumstellar radiative-transfer fits.',
                     'No contemporaneous g-r light curve exists in the retrieved target photometry; optical-infrared changes are an alternative color baseline.',
                     'NEOWISE ends in 2024; it cannot measure the 2025-2026 infrared recovery.'],input_sha256=hashes)
    (out/'dust_analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(4,1,figsize=(12,11),sharex=True,layout='constrained')
    for camera in sorted({r['camera'] for r in optical}):
        data=[r for r in optical if r['camera']==camera]
        base=[r['mag'] for r in data if r['year']==2018]
        if len({r['day'] for r in data if r['year']==2018})<10:continue
        groups=defaultdict(list)
        for r in data:groups[int((r['jd']-2458119.5)//30)].append(r)
        x=[];y=[]
        for g in groups.values():
            if len({r['day'] for r in g})<3:continue
            x.append(float(Time(np.median([r['jd'] for r in g]),format='jd').decimalyear))
            y.append(float(np.median([r['mag'] for r in g])-np.median(base)))
        axes[0].plot(x,y,'.-',lw=.8,ms=3,label=camera,alpha=.7)
    axes[0].set_ylabel('ASAS-SN g change\n(mag, camera offsets)');axes[0].legend(ncol=6,fontsize=8)
    for i,band in enumerate(['W1','W2'],1):
        data=results['primary']['bands'][band]
        for label,style in [('raw','o--'),('corrected','o-')]:
            v=data[label+'_visits'];base=np.median([r['mag'] for r in v if r['year']==2018])
            axes[i].plot([r['decimal_year'] for r in v],[r['mag']-base for r in v],style,ms=4,lw=1,label=label+' visit median')
        axes[i].set_ylabel(band+' change\n(mag, 2018 reference)');axes[i].legend(fontsize=8)
    for label,style in [('raw','o--'),('corrected','o-')]:
        v=results['primary']['colors'][label]['visits'];base=np.median([r['mag'] for r in v if r['year']==2018])
        axes[3].plot([r['decimal_year'] for r in v],[r['mag']-base for r in v],style,ms=4,lw=1,label=label)
    axes[3].set_ylabel('W1 - W2 change\n(mag)');axes[3].set_xlabel('Year');axes[3].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2);ax.axhline(0,color='gray',lw=.6)
        ax.axvspan(2018,2019,alpha=.06,color='green');ax.axvspan(2022,2023,alpha=.06,color='red')
    for ax in axes[:3]:ax.invert_yaxis()
    fig.suptitle('J0754: optical fading versus native NEOWISE visits\nInfrared seasonal/calibration effects require comparison-star correction; no causal fit')
    fig.savefig(out/'optical_infrared.png',dpi=170);plt.close(fig)
    for variant,v in results.items():
        print(variant,{b:(z['target_raw_frames'],z['target_corrected_frames'],z['raw_contrast'],z['corrected_contrast']) for b,z in v['bands'].items()})
    print('models',models)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--asassn',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.previous,a.asassn,a.out)
