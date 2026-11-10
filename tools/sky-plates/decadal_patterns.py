"""Retrospective, passband-separated J0754 event and recurrence checks."""
import csv, hashlib, json
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from asassn_replication import parse_row

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-decadal-patterns-20261109'
PILOT=ROOT/'results/sky-plates-j0754-20261107'
HIST=ROOT/'results/sky-plates-j0754-models-20261108/history'

def read(path):
    return list(csv.DictReader(path.open(newline='')))

def annual(rows):
    groups=defaultdict(list)
    for r in rows: groups[int(float(r['year']))].append(r)
    return [dict(year=y,time=float(np.median([float(r['year']) for r in g])),n=len(g),days=len({r.get('day',str(r['year'])) for r in g}),
                 mag=float(np.median([float(r['corrected_mag']) for r in g])),
                 series=sorted({r['series'] for r in g})) for y,g in sorted(groups.items())]

def contrast(rows, windows, series=None, minimum_annual_n=1):
    if series: rows=[r for r in rows if r['series']==series]
    bins=[r for r in annual(rows) if r['n']>=minimum_annual_n]
    selected=[[r for r in bins if a<=r['year']<b] for a,b in windows]
    if any(len(g)<2 for g in selected):return dict(status='INSUFFICIENT',years=[len(g) for g in selected])
    m=[float(np.mean([r['mag'] for r in g])) for g in selected]
    t=[float(np.mean([r['time'] for r in g])) for g in selected]
    baseline=float(np.interp(t[1],[t[0],t[2]],[m[0],m[2]]))
    return dict(status='MEASURED',annual_equal_weight_mags=m,times=t,years=[len(g) for g in selected],
                points=[sum(r['n'] for r in g) for g in selected],
                rise=m[1]-m[0],recovery=m[2]-m[1],peak_minus_linear_baseline=m[1]-baseline)

def template_fit(rows, template_times, template_mags, shift, trend=False):
    times=np.array([float(r['year']) for r in rows]); use=(times>=min(template_times)+shift)&(times<=max(template_times)+shift)
    selected=[r for r,u in zip(rows,use) if u]; t=times[use]
    y=np.array([float(r['corrected_mag']) for r in selected])
    z=np.interp(t-shift,template_times,template_mags)
    series=sorted({r['series'] for r in selected})
    x=np.array([[float(r['series']==s) for s in series] for r in selected])
    if trend:x=np.column_stack([x,t-np.mean(t)])
    model=np.column_stack([x,z]); b=np.linalg.lstsq(model,y,rcond=None)[0]
    nullb=np.linalg.lstsq(x,y,rcond=None)[0]; oneb=np.linalg.lstsq(x,y-z,rcond=None)[0]
    return dict(n=len(y),time_range=[float(min(t)),float(max(t))],series=series,
                amplitude=float(b[-1]),rms_free=float(np.sqrt(np.mean((y-model@b)**2))),
                rms_no_repeat=float(np.sqrt(np.mean((y-x@nullb)**2))),
                rms_equal_depth_repeat=float(np.sqrt(np.mean((y-z-x@oneb)**2))),
                rank=int(np.linalg.matrix_rank(model)),columns=model.shape[1],condition=float(np.linalg.cond(model)),
                excluded_outside_template=int((~use).sum()))

def main():
    OUT.mkdir(exist_ok=True)
    files=[PILOT/'dasch_differential.csv',HIST/'dasch_all_rows_quality_audit.csv',
           HIST/'controls/annual_exact_frame_differential.csv',
           ROOT/'results/sky-plates-j0754-independent-20261108/asassn_v1_lightcurve.csv']
    hist=read(files[0]); a=annual(hist); raw=read(files[1]); controls=[]
    targetkeys={r['key'] for r in hist}
    windows={'1940s':[(1941,1945),(1945,1949),(1949,1957)],'1980s':[(1969,1976),(1982,1987),(1987,1990)]}
    details={name:dict(windows=w,pooled=contrast(hist,w),by_series={s:contrast(hist,w,s) for s in sorted({r['series'] for r in hist})}) for name,w in windows.items()}
    for name,w in windows.items():
        details[name]['minimum3_sensitivity']={s:contrast(hist,w,s,3) for s in ['ac','rh','dnb']}
    for p in sorted(PILOT.glob('control_*_placebo.csv')):
        files.append(p); cr=[r for r in read(p) if r['key'] in targetkeys]
        controls.append(dict(name=p.stem,features={name:dict(pooled=contrast(cr,w),by_series={s:contrast(cr,w,s) for s in ['ac','rh','dnb']},minimum3_sensitivity={s:contrast(cr,w,s,3) for s in ['ac','rh','dnb']}) for name,w in windows.items()}))
    coverage=[]
    for year in range(1889,1991):
        rr=[r for r in raw if int(r['calendar_year'])==year]
        coverage.append(dict(year=year,rows=len(rr),primary=sum(r['primary_quality_pass']=='True' for r in rr),retained=sum(r['differential_retained']=='True' for r in rr),no_detection=sum('no_photometric_detection' in r['rejection_reasons'] for r in rr)))
    t=sorted(float(r['year']) for r in hist)
    gaps=sorted([dict(start=x,end=y,years=y-x) for x,y in zip(t,t[1:])],key=lambda r:-r['years'])[:12]
    template=[r for r in a if 1932<=r['year']<=1943]
    tt=np.array([r['time'] for r in template]); yy=np.array([r['mag'] for r in template]);yy-=np.median([r['mag'] for r in template if r['year']>=1941])
    late=[r for r in hist if 1969<=float(r['year'])<1990]
    fits={str(shift):{str(trend):template_fit(late,tt,yy,shift,trend) for trend in [False,True]} for shift in [43.,44.5,46.]}
    asas=read(files[2]); modern=[]
    for r in read(files[3]):
        q=parse_row(r)
        if q['accepted']:
            modern.append(dict(year=float(Time(q['jd'],format='jd').decimalyear),corrected_mag=q['mag'],series=q['camera']+' '+q['band'],day=q['day']))
    ma={s:annual([r for r in modern if r['series']==s]) for s in sorted({r['series'] for r in modern})}
    modern_summary={s:dict(annual=aa,faintest_well_sampled=max([r for r in aa if r['days']>=10],key=lambda r:r['mag'])) for s,aa in ma.items() if any(r['days']>=10 for r in aa)}
    anchor_old=next(r['time'] for r in a if r['year']==1933)
    anchor_new=modern_summary['bq g']['faintest_well_sampled']['time']
    proposed=(anchor_new-anchor_old)/6
    schedule=[]
    for cycle in range(-3,7):
        center=anchor_old+cycle*proposed
        hh=[r for r in hist if abs(float(r['year'])-center)<=3]
        center_jd=float(Time(center,format='decimalyear',scale='tt').jd)
        rr=[r for r in raw if abs(float(r['date_jd'])-center_jd)<=3*365.25]
        schedule.append(dict(cycle=cycle,center=center,window=[center-3,center+3],historical_retained=len(hh),historical_catalogue_rows=len(rr),annual=annual(hh)))
    cyclefits={str(cycle):{str(trend):template_fit(hist,tt,yy,cycle*proposed,trend) for trend in [False,True]} for cycle in [1,2,3]}
    result=dict(status='RETROSPECTIVE_DESCRIPTIVE',historical_annual=a,coverage=coverage,largest_retained_gaps=gaps,
                local_features=details,matched_controls=controls,half_interval_template=dict(shift_years=44.5,template_times=tt.tolist(),template_mags=yy.tolist(),fits=fits),
                legacy_asas=asas,modern=modern_summary,exploratory_15year=dict(period=proposed,anchors=[anchor_old,anchor_new],cycles_between_anchors=6,schedule=schedule,template_fits=cyclefits),input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,ax=plt.subplots(3,1,figsize=(12,10),layout='constrained',gridspec_kw={'height_ratios':[2,1,2]})
    baseline=13.114900375
    ax[0].scatter([float(r['year']) for r in hist],[float(r['corrected_mag'])-baseline for r in hist],s=6,alpha=.15,color='navy')
    for sufficient,marker in [(True,'o'),(False,'o')]:
        b=[r for r in a if (r['n']>=3)==sufficient]
        ax[0].scatter([r['time'] for r in b],[r['mag']-baseline for r in b],s=25,facecolors='navy' if sufficient else 'none',edgecolors='navy',label='Annual median, >=3 detections' if sufficient else 'Only 1-2 detections')
    ax[0].set(xlim=(1888,1991),ylabel='DASCH relative photographic magnitude',title='Historical observations: major 1930s fading plus smaller variations')
    ax[0].invert_yaxis();ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
    ax[1].bar([r['year']+.5 for r in coverage],[r['rows'] for r in coverage],color='lightgray',label='Catalogue rows, not all calibrated detections')
    ax[1].bar([r['year']+.5 for r in coverage],[r['retained'] for r in coverage],color='navy',label='Retained differential detections')
    ax[1].set(xlim=(1888,1991),ylabel='Records per year',xlabel='Calendar year');ax[1].legend(fontsize=8)
    ax[2].plot([float(r['year'])+.5 for r in asas],[float(r['target_corrected_median'])-12.153 for r in asas],'s-',label='ASAS V, relative to 2004',color='purple')
    for s,rows in ma.items():
        base=[r['mag'] for r in rows if r['year']==2018 and r['n']>=20]
        if not base:continue
        rows=[r for r in rows if r['n']>=20]
        ax[2].plot([r['time'] for r in rows],[r['mag']-base[0] for r in rows],'.-',label='ASAS-SN '+s+', relative to 2018')
    ax[2].set(xlim=(2001,2027),ylabel='Within-stream magnitude change',xlabel='Calendar year',title='Modern episodes remain in separate filters and camera zero points');ax[2].invert_yaxis();ax[2].grid(alpha=.2);ax[2].legend(fontsize=7,ncol=2)
    fig.savefig(OUT/'decadal_pattern.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(features=details,template=fits,gaps=gaps[:4],modern={s:r['faintest_well_sampled'] for s,r in modern_summary.items()}),indent=2))

if __name__=='__main__':main()
