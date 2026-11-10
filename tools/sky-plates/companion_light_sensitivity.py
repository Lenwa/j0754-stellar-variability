"""Conditional independent-epoch CK04 mixture shape sensitivity; no likelihood."""
import json, hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from compare_spherex import source_filename
from compare_stellar_models import dust_grid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-companion-light-20261109'
MODELS=ROOT/'results/sky-plates-j0754-models-20261108/models'
SPH=ROOT/'results/sky-plates-j0754-spherex-20261108'
FRACTIONS=np.arange(0,.901,.05)
AV=np.arange(0,1.501,.1)
TEMPS=[3500,4000,4500,5000,5500,6000,6500]

def mix(primary,companion,fraction):
    return (1-np.asarray(fraction))*np.asarray(primary)+np.asarray(fraction)*np.asarray(companion)

def shape_score(observed,prediction):
    if np.any(np.asarray(observed)<=0) or np.any(np.asarray(prediction)<=0): raise ValueError('Positive inputs required')
    residual=np.log(np.asarray(observed)/np.asarray(prediction))
    residual=residual-residual.mean(axis=-1,keepdims=True)
    return np.sqrt(np.mean(residual**2,axis=-1))

def datasets(full=False):
    xpfile=ROOT/'results/sky-plates-j0754-models-20261108/spectra/gaia_xp_sampled_calibrated.csv'
    xp=np.loadtxt(xpfile,delimiter=',',skiprows=1)
    centers=np.arange(415,1000,30); obs=[]; sample=[]
    for center in centers:
        mask=(xp[:,0]>=center-15)&(xp[:,0]<center+15)
        obs.append(np.mean(xp[mask,1])); sample.append(xp[mask,0]*10)
    # Endpoint1000 excluded by final bin; equal20 broad bins, no independent error claim.
    data=[dict(name='Gaia_XP_mean',x=np.asarray(sample),obs=np.array(obs),kind='flambda',epoch='DR3 multi-epoch mean',source=str(xpfile))]
    lookup={r['file']:r for r in json.loads((SPH/'local/local_extraction.json').read_text())['records']}
    rows=[]
    for n in Table.read(SPH/'native_IRSA_photometry.tbl',format='ascii.ipac'):
        r=lookup.get(source_filename(n['obs_publisher_did']))
        if r and r['primary_quality'] and int(n['flags'])==1<<21 and float(n['fit_ql'])<5:
            rows.append(dict(r,nwave=float(n['lambda']),nwidth=float(n['lambda_width']),nflux=float(n['flux'])))
    for red in ['native','local']:
        for season in ['2025','2026']:
            rr=[r for r in rows if r['season']==season]
            w=np.array([r['nwave'] if red=='native' else r['wavelength_um'] for r in rr])
            width=np.array([r['nwidth'] if red=='native' else r['width_um'] for r in rr])
            obs=np.array([r['nflux'] if red=='native' else r['flux_ujy'] for r in rr])
            sel=(w>=.75)&(w<=(5.0 if full else 2.4))&(obs>0)
            data.append(dict(name=f'SPHEREx_{red}_{season}',x=(w[sel,None]+width[sel,None]*np.linspace(-.5,.5,31))*1e4,obs=obs[sel],kind='fnu',epoch=season,source=str(SPH)))
    return data

def run(wide=False):
    grain=dust_grid(ROOT/'results/sky-plates-j0754-dust-20261108/draine_rv31_extinction.txt')
    models={}; manifest=[]
    for path in sorted(MODELS.glob('ck*.fits')):
        metal,temp=path.stem.split('_'); temp=int(temp)
        with fits.open(path) as h:
            wave=np.asarray(h[1].data['WAVELENGTH'],float)
            for g in ['g15','g20','g25','g30','g45']:
                if g in h[1].columns.names:
                    flux=np.asarray(h[1].data[g],float)
                    if np.all(flux[wave>=3500]>=0) and np.interp(1e4,wave,flux)>0:
                        models[(metal,temp,g)]=(wave,flux/np.interp(1e4,wave,flux))
        manifest.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    results=[]; fig,axs=plt.subplots(3,2,figsize=(12,12),layout='constrained')
    for data,ax in zip(datasets(),axs.flat):
        x=data['x']; obs=data['obs']; cache={}
        opacity=np.exp(np.interp(np.log(x/1e4),np.log(grain[:,0]),np.log(grain[:,3]))-np.interp(np.log(.55),np.log(grain[:,0]),np.log(grain[:,3])))
        extinction=10**(-.4*AV[:,None,None]*opacity)
        for key,(w,f) in models.items():
            vals=np.interp(x,w,f)
            if data['kind']=='fnu':vals=vals*x*x/2.99792458e18
            cache[key]=np.mean(vals[None,:,:]*extinction,axis=-1)
        profile=np.full((len(TEMPS),len(FRACTIONS)),np.inf)
        best_single=None; best_mix=None; records=[]
        for key,pred in cache.items():
            metal,temp,g=key
            if g=='g45' or not (3500 if wide else 4000)<=temp<=(6500 if wide else 5500):continue
            sc=shape_score(obs,pred); ia=int(np.argmin(sc))
            if best_single is None or sc[ia]<best_single['rms']:
                best_single=dict(model=list(key),av=float(AV[ia]),rms=float(sc[ia])); single_curve=pred[ia]*np.exp(np.mean(np.log(obs/pred[ia])))
            for it,t2 in enumerate(TEMPS):
                key2=(metal,t2,'g45')
                if key2 not in cache:continue
                comp=cache[key2]
                pp=(1-FRACTIONS[:,None,None])*pred[None,:,:]+FRACTIONS[:,None,None]*comp[None,:,:]
                score=shape_score(obs,pp); minav=score.min(axis=1)
                profile[it]=np.minimum(profile[it],minav)
                jf,ja=np.unravel_index(np.argmin(score),score.shape)
                if best_mix is None or score[jf,ja]<best_mix['rms']:
                    best_mix=dict(primary=list(key),companion=list(key2),fraction_1um=float(FRACTIONS[jf]),av=float(AV[ja]),rms=float(score[jf,ja])); mix_curve=pp[jf,ja]*np.exp(np.mean(np.log(obs/pp[jf,ja])))
        envelopes={}
        for tol in [.03,.05,.10]:
            envelopes[str(tol)]={str(t):([float(FRACTIONS[j]) for j in np.flatnonzero(profile[i]<=tol)]) for i,t in enumerate(TEMPS)}
        summary=dict(name=data['name'],epoch=data['epoch'],n=len(obs),single=best_single,mixture=best_mix,temperatures=TEMPS,fractions=FRACTIONS.tolist(),profile_rms=profile.tolist(),compatible_fractions=envelopes,identical_spectrum_degenerate=True)
        results.append(summary)
        center=x.mean(axis=1)/1e4; order=np.argsort(center)
        ax.plot(center,obs/np.median(obs),'.',label='observed')
        ax.plot(center[order],(single_curve/np.median(obs))[order],label='giant')
        ax.plot(center[order],(mix_curve/np.median(obs))[order],label='mixture')
        ax.set_title(data['name']);ax.set_xlabel('Wavelength (micron)');ax.set_ylabel('Scaled flux');ax.legend(fontsize=8)
        print(data['name'],best_single,best_mix,flush=True)
    axs.flat[-1].axis('off');fig.savefig(OUT/('shape_comparison_wide.png' if wide else 'shape_comparison.png'),dpi=150)
    (OUT/('results_wide.json' if wide else 'results.json')).write_text(json.dumps(results,indent=2,allow_nan=False))
    for path in [ROOT/'results/sky-plates-j0754-models-20261108/spectra/gaia_xp_sampled_calibrated.csv',SPH/'local/local_extraction.json',SPH/'native_IRSA_photometry.tbl',ROOT/'results/sky-plates-j0754-dust-20261108/draine_rv31_extinction.txt']:
        manifest.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (OUT/'input_hashes.json').write_text(json.dumps(manifest,indent=2))
    longwave=[]
    for data, result in zip(datasets(full=True),results):
        if data['kind']!='fnu':continue
        x=data['x']; centers=x.mean(axis=1)/1e4;obs=data['obs'];train=centers<=2.4
        opacity=np.exp(np.interp(np.log(x/1e4),np.log(grain[:,0]),np.log(grain[:,3]))-np.interp(np.log(.55),np.log(grain[:,0]),np.log(grain[:,3])))
        for which in ['single','mixture']:
            b=result[which]; key=tuple(b['model'] if which=='single' else b['primary']);w,f=models[key]
            flux=np.interp(x,w,f)
            if which=='mixture':
                w2,f2=models[tuple(b['companion'])];frac=b['fraction_1um'];flux=(1-frac)*flux+frac*np.interp(x,w2,f2)
            pred=np.mean(flux*x*x/2.99792458e18*10**(-.4*b['av']*opacity),axis=1)
            scale=np.exp(np.mean(np.log(obs[train]/pred[train])));ratio=obs/(pred*scale)
            bands={}
            for low,high in [(2.4,3.75),(3.8,4.4),(4.4,5.0)]:
                select=(centers>low)&(centers<=high)
                bands[f'{low}_{high}']=dict(n=int(select.sum()),median_observed_model=float(np.median(ratio[select])),log_rms_fixed_scale=float(np.sqrt(np.mean(np.log(ratio[select])**2))))
            longwave.append(dict(dataset=data['name'],model=which,bands=bands))
    (OUT/('longwave_wide.json' if wide else 'longwave.json')).write_text(json.dumps(longwave,indent=2,allow_nan=False))
if __name__=='__main__':
    import sys
    run('--wide' in sys.argv)


