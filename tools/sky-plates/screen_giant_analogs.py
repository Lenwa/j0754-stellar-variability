"""Exploratory giant comparison selection; catalogue rank precedes archive outcomes."""
import csv
import math
import json
import hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy.io.votable import parse_single_table
from analysis_core import quality_dasch, separation_arcsec
from summarize_modern import good_ztf

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/giant-screening-20261111'
RECON=ROOT/'results/sky-plates-recon-20261107'
REFERENCE=dict(teff=5205.7764,luminosity=54.640694)
REFERENCE['radius']=math.sqrt(REFERENCE['luminosity'])*(5772/REFERENCE['teff'])**2


def finite(value):
    try:
        result=float(value)
        return result if math.isfinite(result) else None
    except (TypeError,ValueError):return None


def parse_vizier_tsv(raw):
    """ASU response for exactly one table; retain missing fields and reject errors."""
    if '#INFO\tError=' in raw or 'OVERFLOW' in raw:
        raise ValueError('Archive error or truncated result')
    lines=[line for line in raw.splitlines() if line.strip() and not line.startswith('#')]
    if len(lines)<3 or len(set(lines[2].replace('\t','').replace('-','').strip())):
        raise ValueError('Missing ASU table header, unit or separator row')
    headers=lines[0].split('\t')
    result=[]
    for line in lines[3:]:
        fields=line.split('\t')
        if len(fields)!=len(headers):raise ValueError('Unexpected row width')
        result.append(dict(zip(headers,[v.strip() for v in fields])))
    return result


def analog_score(row):
    t,r,l=[finite(row.get(k)) for k in ('Teff','Rad','Lbol')]
    if any(v is None or v<=0 for v in (t,r,l)):return None
    return ((t-REFERENCE['teff'])/500)**2+(math.log10(r/REFERENCE['radius'])/.25)**2+(math.log10(l/REFERENCE['luminosity'])/.35)**2


def rank_catalog(rows):
    result=[dict(r,physical_score=analog_score(r),physical_rank=None) for r in rows]
    eligible=sorted([r for r in result if r.get('LumClass')=='giant' and r['physical_score'] is not None],key=lambda r:(r['physical_score'],r['KIC']))
    for rank,row in enumerate(eligible,1):row['physical_rank']=rank
    return result


def annual_nights(times,magnitudes,minimum=10):
    """Daily medians then year medians; no repeated-night weighting or interpolation."""
    nights=defaultdict(list)
    for t,m in zip(times,magnitudes):
        if np.isfinite(t) and np.isfinite(m):nights[math.floor(t)].append(float(m))
    years=defaultdict(list)
    for day,values in sorted(nights.items()):
        years[int(Time(day+.5,format='mjd').decimalyear)].append(float(np.median(values)))
    return [dict(year=y,nights=len(values),median_mag=float(np.median(values)),p16=float(np.quantile(values,.16)),p84=float(np.quantile(values,.84))) for y,values in sorted(years.items()) if len(values)>=minimum]


def gap_series(t,y,max_gap):
    t,y=np.asarray(t,float),np.asarray(y,float)
    indices=np.flatnonzero(np.diff(t)>max_gap)+1
    return np.insert(t,indices,np.nan),np.insert(y,indices,np.nan)


def write_csv(path,rows):
    if not rows:return
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)


def rotation_screen(t,y,e):
    """Descriptive joint linear-trend/sinusoid scan, with at least two trial cycles."""
    t,y,e=map(lambda a:np.asarray(a,float),(t,y,e))
    keep=np.isfinite(t)&np.isfinite(y)&np.isfinite(e)&(e>0)
    t,y,e=t[keep],y[keep],e[keep]
    if len(t)<30 or np.ptp(t)<30:return None
    max_period=min(100.,np.ptp(t)/2)
    f=np.linspace(1/max_period,1/5,3601)
    n=np.column_stack([np.ones(len(t)),(t-np.mean(t))/np.ptp(t)])/e[:,None]
    q=np.linalg.qr(n,mode='reduced')[0]
    residual=y/e-q@(q.T@(y/e));denom=residual@residual
    if denom<1e-20:return None
    power=[]
    for frequency in f:
        angle=2*np.pi*(t-t.min())*frequency
        a=np.column_stack([np.sin(angle),np.cos(angle)])/e[:,None]
        a-=q@(q.T@a)
        b=np.linalg.qr(a,mode='reduced')[0]
        power.append(float(np.sum((b.T@residual)**2)/denom))
    index=int(np.argmax(power))
    return dict(best_period_days=float(1/f[index]),power=power[index],period_range_days=[5.,float(max_period)],points=len(t),span_days=float(np.ptp(t)),status='Descriptive peak; no red-noise/alias validation or orbital identification')


def period_components(raw):
    values=[finite(part.strip()) for part in raw.split(',')]
    return [v for v in values if v is not None and v>0] if all(v is not None and v>0 for v in values) else []


def literature_rows():
    tables=json.loads((OUT/'olah2021_tables.json').read_text(encoding='utf8'))['tables']
    main=next(t for t in tables if t[0]==['KIC No.','Kpmag','Teff [K]','R/R\u2299','lum. class','Prot [d]','no. of flares','binarity'])
    oscillating={r[0]:r[-1] for r in tables[0][1:]}
    return {r[0]:dict(published_rotation_days=finite(r[5]),published_rotation_raw=r[5],published_rotation_components_days=period_components(r[5]),published_binary=r[7],oscillation_note=oscillating.get(r[0],''),Kp=finite(r[1])) for r in main[1:] if len(r)==8}


def archive_summary(target):
    name=target['name'];result=dict(target=target)
    path=OUT/(name+'_gaia.json')
    if path.exists():
        data=json.loads(path.read_text(encoding='utf8'));cols=[m['name'] for m in data['metadata']]
        rows=[dict(zip(cols,row)) for row in data['data']]
        for row in rows:row['separation_arcsec']=separation_arcsec(target['ra'],target['dec'],row['ra'],row['dec'])
        if rows:
            nearest=min(rows,key=lambda r:r['separation_arcsec'])
            result['gaia_nearest']=nearest
            result['gaia_neighbors_21arcsec']=sum(r['source_id']!=nearest['source_id'] and r['separation_arcsec']<21 for r in rows)
            result['gaia_field_rows']=len(rows)
            result['identity_note']='Catalogue position match; J2000 input versus Gaia epoch 2016, no resolved component attribution'
    path=OUT/(name+'_mast.json')
    if path.exists():
        data=json.loads(path.read_text(encoding='utf8'));result['mast_status']=data.get('status')
        counts=defaultdict(int)
        for r in data.get('data',[]):counts[str(r.get('obs_collection'))+' / '+str(r.get('dataproduct_type'))]+=1
        result['mast_counts']=dict(counts)
        result['mast_paging']=data.get('paging')
    return result


def tang_curves():
    results={}
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for i,name in enumerate(['J0736','J0830']):
        with (RECON/(name+'_apass_lightcurve.csv')).open() as stream:raw=list(csv.DictReader(stream))
        result=dict(historical_raw=len(raw));ax=axes[i,0]
        for strict in (False,True):
            rows=[r for r in raw if quality_dasch(r,strict=strict)]
            t=np.array([float(r['date_jd'])-2400000.5 for r in rows]);m=np.array([float(r['magcal_magdep']) for r in rows])
            annual=annual_nights(t,m,minimum=3)
            label='strict' if strict else 'standard'
            result['historical_'+label]=dict(points=len(rows),annual=annual)
            if not strict:
                ax.scatter(Time(t,format='mjd').decimalyear,m,s=3,alpha=.18,color='gray')
            ax.plot(*gap_series([r['year']+.5 for r in annual],[r['median_mag'] for r in annual],1),'.-',lw=.8,label=label)
        ax.set(title=name+': DASCH screening curve',ylabel='Photographic magnitude',xlabel='Year');ax.invert_yaxis();ax.legend(fontsize=8)
        table=parse_single_table(RECON/(name+'_ztf.xml')).to_table()
        raw=[{key:str(row[key]) for key in table.colnames} for row in table]
        accepted=[r for r in raw if good_ztf(r)];groups=defaultdict(list)
        for r in accepted:groups[(r['oid'],r['filtercode'],r['field'],r['ccdid'],r['qid'])].append(r)
        modern=[];seasonal=[];ax=axes[i,1]
        for group,rr in sorted(groups.items(),key=lambda item:-len(item[1])):
            t=np.array([float(r['mjd']) for r in rr]);m=np.array([float(r['mag']) for r in rr]);e=np.array([float(r['magerr']) for r in rr])
            annual=annual_nights(t,m);record=dict(stream=list(group),points=len(rr),annual=annual)
            record['annual_range_mag']=float(np.ptp([a['median_mag'] for a in annual])) if len(annual)>=3 else None
            modern.append(record)
            if len(modern)<=2:
                ax.scatter(Time(t,format='mjd').decimalyear,m-np.median(m),s=3,alpha=.18)
                ax.plot(*gap_series([r['year']+.5 for r in annual],[r['median_mag']-np.median(m) for r in annual],1),'.-',label=group[1]+' / '+group[0][-5:])
            # Only the most observed stream in each filter enters the exploratory seasonal scan.
            if any(s['filter']==group[1] for s in seasonal):continue
            dates=Time(t,format='mjd').to_value('iso',subfmt='date')
            seasons=np.array([int(d[:4])-(int(d[5:7])<7) for d in dates])
            records=[]
            for year in sorted(set(seasons)):
                mask=seasons==year
                # Median duplicate observations per UTC night before period screening.
                nights=defaultdict(list)
                for ti,mi,ei in zip(t[mask],m[mask],e[mask]):nights[math.floor(ti)].append((ti,mi,ei))
                bins=np.array([np.median(v,axis=0) for v in nights.values()])
                if len(bins)<30 or np.ptp(bins[:,0])<120:continue
                peak=rotation_screen(bins[:,0],bins[:,1],np.maximum(.02,bins[:,2]))
                if peak:records.append(dict(season=str(year)+'/'+str(year+1)[2:],**peak))
            seasonal.append(dict(filter=group[1],stream=list(group),seasons=records))
        result['modern']=dict(snapshot='Retained 7 November 2026 ZTF query; fresh 11 November requests timed out.',raw_points=len(raw),quality_points=len(accepted),streams=modern,seasonal_peaks=seasonal)
        results[name]=result
        ax.set(title=name+' ZTF: separate object/filter/detector streams',ylabel='Magnitude relative to stream median',xlabel='Year');ax.invert_yaxis()
        if groups:ax.legend(fontsize=8)
        else:
            ax.set_title(name+': modern archive gap');ax.set_axis_off();ax.text(.5,.5,'No ZTF catalogue rows in retained query.\nFresh request timed out.\nModern behavior remains unmeasured here.',ha='center',va='center',transform=ax.transAxes)
    fig.suptitle('Tang comparison stars: screening curves, without matched-control calibration')
    fig.savefig(OUT/'tang_comparison_lightcurves.png',dpi=150);plt.close(fig)
    return results


def kepler_curves(selected):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(len(selected),2,figsize=(12,12),layout='constrained')
    records=[]
    for i,star in enumerate(selected):
        paths=sorted(OUT.glob('kplr'+star['KIC'].zfill(9)+'-*_llc.fits'))
        for j,path in enumerate(paths):
            with fits.open(path) as hdul:
                data=hdul[1].data;header=hdul[1].header;primary=hdul[0].header
                mask=(data['SAP_QUALITY']==0)&np.isfinite(data['TIME'])
                mask &= np.isfinite(data['SAP_FLUX'])&(data['SAP_FLUX']>0)&np.isfinite(data['PDCSAP_FLUX'])&(data['PDCSAP_FLUX']>0)
                t=np.asarray(data['TIME'][mask],float);sap=np.asarray(data['SAP_FLUX'][mask],float);pdc=np.asarray(data['PDCSAP_FLUX'][mask],float)
                sap=sap/np.median(sap);pdc=pdc/np.median(pdc)
                record=dict(kic=star['KIC'],file=path.name,quarter=primary.get('QUARTER'),raw_points=len(data),accepted_points=len(t),span_days=float(np.ptp(t)),crowdsap=header.get('CROWDSAP'),pdc_p05_p95_percent=((np.quantile(pdc,[.05,.95])-1)*100).tolist())
            binned=defaultdict(list)
            for ti,si,pi in zip(t,sap,pdc):binned[math.floor(ti/.2)].append((ti,si,pi))
            bins=np.array([np.median(v,axis=0) for v in binned.values()])
            record['sap_rotation']=rotation_screen(bins[:,0],bins[:,1],np.ones(len(bins)))
            record['pdc_rotation']=rotation_screen(bins[:,0],bins[:,2],np.ones(len(bins)))
            write_csv(OUT/(path.stem+'_screen.csv'),[dict(bkjd=float(ti),sap_relative=float(si),pdc_relative=float(pi)) for ti,si,pi in zip(t,sap,pdc)])
            records.append(record)
            if j<2:
                ax=axes[i,j];ax.scatter(t-t.min(),(pdc-1)*100,s=1,alpha=.6,color='navy');ax.plot(*gap_series(bins[:,0]-t.min(),(bins[:,1]-1)*100,1),lw=.6,color='darkorange',label='SAP, 0.2-day medians')
                ax.set(title=f"KIC {star['KIC']} / Q{record['quarter']} / {len(t)} retained",ylabel='Relative flux (%)',xlabel='Days from quarter start');ax.grid(alpha=.2)
                if i==0:ax.legend(fontsize=7)
    fig.suptitle('Early and later Kepler quarters for the five physical analogues\nBlue: PDC; orange: SAP. Quarters normalized independently; short Q17 also retained in the table.')
    fig.savefig(OUT/'kepler_screening_lightcurves.png',dpi=150);plt.close(fig)
    return records


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows,selected=shortlist();literature=literature_rows()
    for row in rows:
        row.update(literature.get(row['KIC'],{}))
        row['logg_label_tension']=bool(row.get('LumClass')=='giant' and finite(row.get('logg')) is not None and float(row['logg'])>3.5)
    write_csv(OUT/'screened_catalog.csv',rows)
    selected=sorted([r for r in rows if r['physical_rank'] is not None],key=lambda r:r['physical_rank'])[:5]
    write_csv(OUT/'physical_shortlist.csv',selected)
    archive=[archive_summary(t) for t in json.loads((OUT/'targets.json').read_text())]
    orbits=json.loads((OUT/'tang_gaia_orbits.json').read_text());columns=[r['name'] for r in orbits['metadata']]
    orbits=[dict(zip(columns,row)) for row in orbits['data']]
    orbit_keys=['source_id','nss_solution_type','period','period_error','eccentricity','eccentricity_error','semi_amplitude_primary','rv_n_good_obs_primary','significance','goodness_of_fit']
    orbits=[{k:r.get(k) for k in orbit_keys} for r in orbits]
    historical=tang_curves();kepler=kepler_curves(selected)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    valid=[r for r in rows if r.get('LumClass')=='giant' and r['physical_score'] is not None]
    ax.scatter([float(r['Teff']) for r in valid],[float(r['Rad']) for r in valid],color='lightgray',label='Published giant label')
    for r in selected:
        x,y=float(r['Teff']),float(r['Rad']);ax.scatter(x,y,s=35);ax.annotate(' '+r['KIC'],(x,y),fontsize=8)
    ax.scatter(REFERENCE['teff'],REFERENCE['radius'],marker='*',s=160,color='black',label='J0754 adopted model')
    ax.set(xlabel='Effective temperature (K)',ylabel='Radius (solar radii)',title='Physical comparison screen; catalogue values, not uniform new stellar fits');ax.invert_xaxis();ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(OUT/'physical_analogues.png',dpi=160);plt.close(fig)
    receipts=[json.loads(p.read_text()) for p in sorted(OUT.glob('*.receipt.json'))]
    result=dict(reference=REFERENCE,classification_counts={name:sum(r['LumClass']==name for r in rows) for name in ['giant','dwarf','']},rankable_giants=len(valid),logg_label_tensions=[r['KIC'] for r in rows if r['logg_label_tension']],selected=selected,archives=archive,tang_gaia_candidate_orbits=orbits,tang=historical,kepler=kepler,retrieval_states={r['file']:r['state'] for r in receipts},limitations=['Exploratory comparison shortlist; no incidence or novelty claim.','Dixon 7286 rows were not retrieved and are not screened.','Published luminosity labels retained alongside changed stellar parameters.','Kepler period maxima are descriptive, with no red-noise or alias significance.','No new flare search or resolved host assignment on the comparison stars.','DASCH and ZTF curves lack the matched controls used for J0754.','Quarter normalization prevents measurement of absolute multi-year fading.'])
    (OUT/'screening.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')
    hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.suffix in ['.fits','.tsv','.xml']}
    hashes.update({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('fetch_giant_screening.py')]})
    (OUT/'input_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(classification_counts=result['classification_counts'],rankable_giants=len(valid),kepler_files=len(kepler),tang_points={name:(r['historical_standard']['points'],r['modern']['quality_points']) for name,r in historical.items()})))


def shortlist():
    rows=rank_catalog(parse_vizier_tsv((OUT/'olah2022_stars.tsv').read_text(encoding='utf8')))
    selected=sorted([r for r in rows if r['physical_rank'] is not None],key=lambda r:r['physical_rank'])[:5]
    return rows,selected


if __name__=='__main__':main()
