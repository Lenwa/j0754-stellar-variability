"""Offline GALEX response/background sensitivity and official X-ray counterpart check."""
import hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-galex-erosita-20261111'
OLD=ROOT/'results/sky-plates-j0754-multiwavelength-20261110'
TARGET=SkyCoord(118.69117410685,16.694621231895,unit='deg')

def rows(name):
    data=json.loads((OUT/name).read_text())
    if data.get('status')!='COMPLETE':raise ValueError(data.get('msg'))
    return np.array(data['data']['Tables'][0]['Rows'],float)

def flat_response(xi,eta,t,flat,pixsz):
    # gPhoton curvetools.xieta2colrow/hashresponse: first index is xi.
    ix=(np.asarray(xi)/36000/pixsz+flat.shape[0]/2).astype(int)
    iy=(np.asarray(eta)/36000/pixsz+flat.shape[1]/2).astype(int)
    if np.any((ix<0)|(iy<0)|(ix>=flat.shape[0])|(iy>=flat.shape[1])):raise ValueError('Outside flat')
    scale=1.9946352-1.9679445e-9*t+9.3025231e-19*t*t
    scale=np.where(t>=881881215.995,scale*1.018,scale)
    response=flat[ix,iy]*scale
    if np.any(response<=0):raise ValueError('Nonpositive detector response')
    return response

def bin_edges(start,stop,width):
    edges=np.r_[np.arange(start,stop,width),stop]
    if len(edges)>2 and edges[-1]-edges[-2]<width/2:edges=np.delete(edges,-2)
    return edges

def live_time(a,b,ticks):
    # Global event presence in native 0.05s intervals implements shutter coverage.
    # gPhoton dbasetools.empirical_deadtime: dead fraction = -slope*gcr/intercept.
    x=ticks[:,0]*.05
    overlap=np.maximum(0,np.minimum(x+.05,b)-np.maximum(x,a))
    open_time=overlap.sum()
    if open_time<=0:raise ValueError('No open detector time')
    global_counts=np.sum(ticks[:,1]*overlap/.05)
    dead=.000434730599193*(global_counts/open_time)/77.217817988
    if not 0<=dead<1:raise ValueError('Invalid deadtime')
    return open_time*(1-dead),dead,(b-a)-open_time

def constant_test(n,nb,es,eb,alpha,simulations=100000):
    n,nb,es,eb=map(lambda x:np.asarray(x,float),(n,nb,es,eb))
    bg=nb.sum()/eb.sum();rate=n/es-alpha*nb/eb
    source=max(0.,n.sum()/es.sum()-alpha*bg)
    variance=(source+alpha*bg)/es+alpha*alpha*bg/eb
    weights=1/variance;mean=np.sum(rate*weights)/weights.sum()
    statistic=float(np.sum((rate-mean)**2*weights))
    rng=np.random.default_rng(751416)
    a=rng.poisson((source+alpha*bg)*es,(simulations,len(n)))
    b=rng.poisson(bg*eb,(simulations,len(n)))
    simulated=a/es-alpha*b/eb
    m=np.sum(simulated*weights,axis=1)/weights.sum()
    reference=np.sum((simulated-m[:,None])**2*weights,axis=1)
    return dict(rate=rate.tolist(),sigma=np.sqrt(variance).tolist(),mean_rate=float(mean),statistic=statistic,
                null_simulation_p=float((1+np.count_nonzero(reference>=statistic))/(simulations+1)),
                note='Approximate Poisson reference conditional on measured response tracks; source and annulus simulated under constant rates. No trials correction.')

def galex():
    photons=rows('photons_30arcsec.json');t=photons[:,0]/1000
    radius=TARGET.separation(SkyCoord(photons[:,1],photons[:,2],unit='deg')).arcsec
    with fits.open(OUT/'NUV_flat.fits') as h:response=flat_response(photons[:,3],photons[:,4],t,h[0].data,h[0].header['CDELT2'])
    records=[];fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for aperture in [6,9,12]:
        for width in [20,30,40]:
            bins=[]
            for visit,(start,stop) in enumerate([(822929965.495,822930062.495),(822935895.495,822935977.495)]):
                ticks=rows(f'global_ticks_{visit}.json');edges=bin_edges(start,stop,width)
                for a,b in zip(edges[:-1],edges[1:]):
                    time=(t>=a)&(t<b);src=time&(radius<=aperture);bg=time&(radius>=20)&(radius<30)
                    live,dead,shutter=live_time(a,b,ticks)
                    # n / sum(1/response) is the effective response for weighted photons.
                    rs=src.sum()/np.sum(1/response[src]);rb=bg.sum()/np.sum(1/response[bg])
                    bins.append(dict(visit=visit,time_s=float((a+b)/2-start),duration_s=float(b-a),counts=int(src.sum()),background_counts=int(bg.sum()),live_s=float(live),dead_fraction=float(dead),shutter_s=float(shutter),response=float(rs),background_response=float(rb),source_effective_s=float(live*rs),background_effective_s=float(live*rb)))
            alpha=aperture**2/(30**2-20**2)
            result=constant_test([x['counts'] for x in bins],[x['background_counts'] for x in bins],[x['source_effective_s'] for x in bins],[x['background_effective_s'] for x in bins],alpha)
            result.update(aperture_arcsec=aperture,bin_width_s=width,bins=bins);records.append(result)
            if width==20:
                for visit in [0,1]:
                    mask=np.array([x['visit']==visit for x in bins]);x=np.array([b['time_s'] for b in bins])[mask]
                    axes[0,visit].errorbar(x,np.array(result['rate'])[mask],np.array(result['sigma'])[mask],fmt='o-',label=f'{aperture}" aperture',alpha=.75)
                    if aperture==6:
                        axes[1,visit].plot(x,np.array([b['response'] for b in bins])[mask],'-o',label='NUV response')
                        axes[1,visit].plot(x,1-np.array([b['dead_fraction'] for b in bins])[mask],'-s',label='Live fraction')
    for v in [0,1]:
        axes[0,v].set(title=f'Visit {v+1}: 2 February 2006',ylabel='Background-subtracted corrected ct/s');axes[0,v].legend()
        axes[1,v].set(xlabel='Seconds from visit start',ylabel='Detector factor');axes[1,v].legend()
    fig.suptitle('J0754 GALEX NUV | 20-second bins, short tails merged')
    fig.savefig(OUT/'galex_corrected.png',dpi=160);plt.close(fig)
    return dict(photons=int(len(t)),response_range=[float(response.min()),float(response.max())],tests=records,
                limitation='Apertures tested for time-dependent losses; no infinite-aperture correction or absolute UV luminosity. Sparse 179s coverage cannot constrain multi-hour flare rates. Photon-weight-derived flat tracks make null probabilities approximate.')

def erosita():
    t=Table.read(OUT/'erosita_gaia_counterparts.ecsv')
    keys=['IAUNAME','UID','GDR3_source_id','GDR3_RA','GDR3_DEC','GDR3_Xray_proba','NWAY_Separation_GDR3_ERO','NWAY_p_any','NWAY_p_i','NWAY_match_flag','NWAY_purity6','NWAY_completeness6','NWAY_threshold6']
    matches=[]
    for r in t:
        matches.append({k:(r[k].item() if hasattr(r[k],'item') else str(r[k])) for k in keys})
        matches[-1]={k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in matches[-1].items()}
    return dict(counterparts=matches,source='Official eROSITA DR2 Gaia DR3/NWAY catalogue, all rows searched by UID.',
                limitation='Gaia-catalogue counterpart probabilities are conditional on that input catalogue and its priors; an uncatalogued faint emitter or blend is not excluded. DR2 is catalogue-only: no DR2 image or time-resolved flare test.')

def main():
    result=dict(galex=galex(),erosita=erosita())
    result['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.suffix in ['.json','.fits','.ecsv','.txt'] and p.name not in ['results.json','report.txt']}
    result['source_sha256']={name:hashlib.sha256((ROOT/'tools/sky-plates'/name).read_bytes()).hexdigest() for name in ['fetch_galex_erosita_followup.py','analyze_galex_erosita_followup.py']}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(json.dumps(dict(galex=[{k:r[k] for k in ['aperture_arcsec','bin_width_s','mean_rate','null_simulation_p']} for r in result['galex']['tests']],erosita=result['erosita']),indent=2))

if __name__=='__main__':main()
