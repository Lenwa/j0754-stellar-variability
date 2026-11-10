"""Existing-data, band-integrated broader dust sensitivities; no model probabilities."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import brentq
from fetch_stellar_grids import FILTERS
from stellar_models import read_bandpass,band_flux
from compare_stellar_models import dust_grid,extinction_deficits

BANDS=['g','W1','W2','V','r','I','J','Ks']

def cover_ratio(deficit,c0,c1):
    return (1-c1*np.asarray(deficit))/(1-c0*np.asarray(deficit))

def final_cover(deficit_g,c0,loss_g):
    return c0+loss_g*(1-c0*np.asarray(deficit_g))/np.maximum(deficit_g,1e-100)

def fixed_cover_from_columns(d0,d1,ratio):
    return (1-ratio)/(np.asarray(d1)-ratio*d0)

def gray_selective(gray,selective,ratios):
    return gray+selective*np.asarray(ratios)

def select_rows(pred,obs,tol=.02):
    return np.all(np.abs(np.asarray(pred)[:,:2]-obs[:2])<=tol+1e-12,axis=1)

def mags(ratio):return -2.5*np.log10(np.maximum(ratio,1e-250))

def main(root,out):
    out.mkdir(exist_ok=True,parents=True)
    models=root/'results/sky-plates-j0754-models-20261108/models'
    dust=root/'results/sky-plates-j0754-dust-20261108'
    bands=[read_bandpass(models/FILTERS[b][0],FILTERS[b][1]) for b in BANDS]
    grains={x:dust_grid(dust/f'draine_{x}_extinction.txt') for x in ['rv31','rv4','rv5.5']}
    prior=json.loads((root/'results/sky-plates-j0754-models-20261108/model_comparison.json').read_text())
    tau=np.r_[0,np.geomspace(.0001,300,1200)]
    records=[];envelopes={};trajectory=[]
    def save(family,pred,params,obs,meta):
        p=np.asarray(pred)
        if len(p)==0:return
        ok=select_rows(p,obs)
        idx=int(np.argmin(np.sum((p[:,:2]-obs[:2])**2,axis=1)))
        records.append(dict(family=family,**meta,parameters=params[idx],prediction=dict(zip(BANDS,p[idx].tolist())),
                            training_residual=(p[idx,:2]-obs[:2]).tolist(),W2_check_residual=float(p[idx,2]-obs[2])))
        key=meta['season']+':'+family
        e=envelopes.setdefault(key,dict(count=0,bands={},parameters=[]))
        if np.any(ok):
            e['count']+=int(ok.sum())
            for j,b in enumerate(BANDS):
                v=[float(p[ok,j].min()),float(p[ok,j].max())]
                if b in e['bands']:v=[min(v[0],e['bands'][b][0]),max(v[1],e['bands'][b][1])]
                e['bands'][b]=v
            # Full sampled compatible parameters and predictions retained for audit.
            e['parameters'].extend([{**meta,**params[i],"prediction":p[i].tolist()} for i in np.where(ok)[0]])
    for temp in [4500,5000,5500]:
        with fits.open(models/f'ckm05_{temp}.fits') as h:
            wave=np.array(h[1].data['WAVELENGTH'],float);hot=np.array(h[1].data['g25'],float)
        flux=np.array([band_flux(wave,hot,b) for b in bands])
        screens={law:extinction_deficits(wave,hot,bands,grain,tau) for law,grain in grains.items()}
        if temp==5000:
            old=json.loads((root/'results/sky-plates-j0754-chromaticity-20261109/results.json').read_text())
            for camera in old['trajectory']:
                anchor=np.array([camera['training_g'][1]-camera['training_g'][0],camera['training_W1'][1]-camera['training_W1'][0]])
                rr=10**(-.4*anchor)
                for law,d in screens.items():
                    for t0 in [0,.1,.3,1,3,10]:
                        d0=np.array([np.interp(t0,tau,d[:,j]) for j in range(len(BANDS))])
                        denominator=d[:,0]-rr[0]*d0[0]
                        cv=np.divide(1-rr[0],denominator,out=np.full_like(denominator,np.inf),where=denominator>0)
                        valid=(tau>t0)&(cv>0)&(cv<=1)
                        if not valid.any():continue
                        pred=mags((1-cv[valid,None]*d[valid])/(1-cv[valid,None]*d0))
                        j=int(np.argmin(abs(pred[:,1]-anchor[1])))
                        c=cv[valid][j];training_error=float(pred[j,1]-anchor[1])
                        if abs(training_error)>.002:continue
                        # Hold c and initial column fixed, infer each visit column from g alone.
                        gcurve=mags((1-c*d[:,0])/(1-c*d0[0]))
                        wcurve=mags((1-c*d[:,1])/(1-c*d0[1]))
                        rows=[];outside=[]
                        for row in camera['rows']:
                            dg=row['g']-camera['training_g'][0]
                            if not gcurve[0]<=dg<=gcurve[-1]:
                                outside.append(dict(date=row['date'],g_change=dg,allowed_g_change=[float(gcurve[0]),float(gcurve[-1])]))
                                continue
                            dw=float(np.interp(dg,gcurve,wcurve))
                            prediction=camera['training_W1'][0]+dw
                            rows.append(dict(date=row['date'],W1_pred=prediction,residual=row['W1']-prediction))
                        trajectory.append(dict(camera=camera['camera'],window=camera['window_days'],law=law,tau0=t0,cover=float(c),tau1=float(tau[valid][j]),training_error=training_error,
                            bright_g_transmission=float(1-c*d0[0]),n=len(rows),n_available=len(camera['rows']),outside_domain=outside,rms=float(np.sqrt(np.mean([r['residual']**2 for r in rows]))),rows=rows))
        for season,obs in [('March',np.array([.312,.1175,.08875])),('October',np.array([.3225,.115,.09]))]:
            ratio=10**(-.4*obs);loss=1-ratio
            for law,d in screens.items():
                meta=dict(T=temp,law=law,season=season)
                for c0 in [0,.1,.3,.5,.7]:
                    c1=final_cover(d[:,0],c0,loss[0])
                    valid=(c1<=1)&(c1>=c0)&(d[:,0]>0)
                    p=mags(cover_ratio(d[valid],c0,c1[valid,None]))
                    save('fixed_tau_variable_cover',p,[dict(c0=c0,c1=float(c),tau=float(t)) for c,t in zip(c1[valid],tau[valid])],obs,meta)
                # Fixed coverage, changing column. Initial tau is nuisance, final tau set by g.
                p=[];pars=[]
                for t0 in [0,.1,.3,1,3,10]:
                    d0=np.array([np.interp(t0,tau,d[:,j]) for j in range(len(BANDS))])
                    denom=d[:,0]-ratio[0]*d0[0]
                    c=np.divide(loss[0],denom,out=np.full_like(denom,np.inf),where=denom>0)
                    valid=(tau>t0)&(c<=1)&(c>0)
                    p.extend(mags((1-c[valid,None]*d[valid])/(1-c[valid,None]*d0)))
                    pars.extend([dict(cover=float(cc),tau0=t0,tau1=float(t1),bright_transmission_g=float(1-cc*d0[0]),bright_transmission_W1=float(1-cc*d0[1])) for cc,t1 in zip(c[valid],tau[valid])])
                save('fixed_cover_variable_tau',p,pars,obs,meta)
                # Gray opacity plus standard selective screen: solve training contrast for tau.
                screen_mag=mags(1-d)
                delta=screen_mag[:,0]-screen_mag[:,1]
                t=brentq(lambda v:np.interp(v,tau,delta)-(obs[0]-obs[1]),0,3)
                spec=np.array([np.interp(t,tau,screen_mag[:,j]) for j in range(len(BANDS))])
                gray=float(obs[0]-spec[0])
                if gray>=0:save('gray_plus_standard_screen',[spec+gray],[dict(tau=t,gray_mag=gray)],obs,meta)
            # Smooth power laws are phenomenological; not a grain-size inference.
            if temp==5000:
                for c0 in [0,.3]:
                    allpred=[];pars=[]
                    for alpha in np.linspace(0,2.5,101):
                        # Reuse integrator with synthetic opacity table.
                        w=np.geomspace(.1,100,3000);g=np.zeros((len(w),6));g[:,0]=w;g[:,3]=(w/.48)**(-alpha)
                        d=extinction_deficits(wave,hot,bands,g,tau)
                        c1=final_cover(d[:,0],c0,loss[0])
                        valid=(c1<=1)&(c1>=c0)&(d[:,0]>0)
                        p=mags(cover_ratio(d[valid],c0,c1[valid,None]))
                        # Keep all training-compatible and per-alpha best; no W2 selection.
                        keep=select_rows(p,obs);keep[np.argmin(np.sum((p[:,:2]-obs[:2])**2,axis=1))]=True
                        allpred.extend(p[keep]);pars.extend([dict(alpha=float(alpha),c0=c0,c1=float(c),tau=float(t)) for c,t in zip(c1[valid][keep],tau[valid][keep])])
                    save('powerlaw_partial_cover',allpred,pars,obs,dict(T=temp,season=season,c0=c0,law='phenomenological'))
                def fullpred(alpha):
                    w=np.geomspace(.1,100,3000);g=np.zeros((len(w),6));g[:,0]=w;g[:,3]=(w/.48)**(-alpha)
                    d=extinction_deficits(wave,hot,bands,g,tau)
                    sm=mags(1-d);t=np.interp(obs[0],sm[:,0],tau)
                    return np.array([np.interp(t,tau,sm[:,j]) for j in range(len(BANDS))]),t
                a=brentq(lambda a:fullpred(a)[0][1]-obs[1],0,2.5)
                pred,t=fullpred(a)
                save('powerlaw_uniform_screen',[pred],[dict(alpha=a,tau=t)],obs,dict(T=temp,season=season,law='phenomenological'))
            # Independent increasing spots plus uniform standard foreground column.
            tsmall=np.r_[0,np.linspace(.001,.5,300)]
            dh=extinction_deficits(wave,hot,bands,grains['rv31'],tsmall)
            htrans=flux*(1-dh)
            for cooltemp in range(3500,temp,250):
                with fits.open(models/f'ckm05_{cooltemp}.fits') as h:cool=np.array(h[1].data['g25'],float)
                cf=np.array([band_flux(wave,cool,b) for b in bands])
                dc=extinction_deficits(wave,cool,bands,grains['rv31'],tsmall);ctrans=cf*(1-dc)
                for f0 in [0,.1,.3]:
                    bright=(1-f0)*flux+f0*cf
                    f1=(htrans[:,0]-ratio[0]*bright[0])/(htrans[:,0]-ctrans[:,0])
                    valid=(f1>=f0)&(f1<=1)
                    faint=(1-f1[valid,None])*htrans[valid]+f1[valid,None]*ctrans[valid]
                    pred=mags(faint/bright)
                    save('spots_plus_standard_screen',pred,[dict(Tspot=cooltemp,f0=f0,f1=float(f),tau=float(t)) for f,t in zip(f1[valid],tsmall[valid])],obs,dict(T=temp,season=season,law='rv31'))
    result=dict(bands=BANDS,observed={'March':[.312,.1175,.08875],'October':[.3225,.115,.09]},records=records,envelopes=envelopes,variable_column_trajectories=trajectory,
        prior_summary=prior['summary'],prior_envelopes=prior['compatibility_envelopes'],
        qualification='All fitting/selection uses g,W1 only. W2 check retrospective, not new validation. Sampled W1+/-0.02mag regions conditional on exact g anchoring, not full two-band tolerance volumes or confidence intervals. Fixed atmosphere nuisance grid; no thermal/scattered dust light.',
        hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*models.glob('*.fits'),*models.glob('*.xml'),*dust.glob('draine*txt'),root/'results/sky-plates-j0754-chromaticity-20261109/results.json',root/'results/sky-plates-j0754-models-20261108/model_comparison.json']})
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    for k,e in envelopes.items():print(k,e['count'],e['bands'].get('W2'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.root,a.out)


