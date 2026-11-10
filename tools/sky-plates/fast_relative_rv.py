"""Bounded relative FAST line-shift diagnostics, not a calibrated orbital fit."""
import json, hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar, least_squares
C=299792.458
ROOT=Path(__file__).resolve().parents[2]

def normalized(w,f,width=70):
    smooth=gaussian_filter1d(f,width/np.median(np.diff(w)),mode='nearest')
    return f/np.maximum(smooth,np.median(smooth)*.05)-1

def shift_fit(w,y,tw,ty,lo,hi):
    keep=(w>=lo)&(w<=hi)&np.isfinite(y)
    x=w[keep];v=y[keep]
    def score(dv):
        model=np.interp(x/(1+dv/C),tw,ty)
        a=np.column_stack([model,np.ones(len(x))])
        b=np.linalg.lstsq(a,v,rcond=None)[0]
        return np.mean((v-a@b)**2)
    grid=np.arange(-150,150.01,5)
    j=int(np.argmin([score(q) for q in grid]))
    best=minimize_scalar(score,bounds=(max(-150,grid[j]-6),min(150,grid[j]+6)),method='bounded')
    model=np.interp(x/(1+best.x/C),tw,ty)
    nuisance=np.linalg.lstsq(np.column_stack([model,np.ones(len(x))]),v,rcond=None)[0]
    return dict(shift_kms=float(best.x),rms=float(np.sqrt(best.fun)),boundary=bool(abs(best.x)>149),n=len(x),template_amplitude=float(nuisance[0]),offset=float(nuisance[1]),positive_template=bool(nuisance[0]>0))

def sky_center(w,y):
    k=(w>5562)&(w<5592); x=w[k]; z=y[k]
    if len(x)<10:return None
    base=np.median(z); amp=float(np.max(z)-base)
    def fun(p):return p[0]+p[1]*(x-5577)+p[2]*np.exp(-.5*((x-p[3])/p[4])**2)-z
    fit=least_squares(fun,[base,0,max(amp,1e-5),5577,2],bounds=([-np.inf,-np.inf,0,5570,.3],[np.inf,np.inf,np.inf,5584,8]))
    p=fit.x; rms=np.sqrt(np.mean(fun(p)**2))
    return dict(center_A=float(p[3]),sigma_A=float(p[4]),amplitude=float(p[2]),residual_rms=float(rms),amplitude_over_residual=float(p[2]/max(rms,1e-20)),boundary=bool(p[3]<5570.01 or p[3]>5583.99 or p[4]>7.99))

def main():
    inp=ROOT/'results/sky-plates-j0754-spherex-20261108/literature/fast'
    out=ROOT/'results/sky-plates-j0754-fast-rv-20261109';out.mkdir(exist_ok=True)
    rows=[]
    for path in sorted(inp.glob('*.fits')):
        with fits.open(path) as hd:
            h=hd[0].header; a=np.array(hd[0].data,float)[:,0,:]
            w=h['CRVAL1']+(np.arange(a.shape[1])+1-h['CRPIX1'])*h['CD1_1']
            rows.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),date=h['DATE-OBS'],mjd=float(h['MJD']),grating=int(h['DISPERSE']),bcv=float(h['BCV']),pipeline_velocity=float(h['VELOCITY']),pipeline_czerr=float(h['CZERR']),pipeline_xcorr_error=float(h['CZXCERR']),w=w,y=normalized(w,a[0]),flux=a[0],sky=sky_center(w,a[2])))
    templates={300:next(r for r in rows if r['file'].startswith('0135.')),600:next(r for r in rows if r['file'].startswith('0089.'))}
    records=[]
    for r in rows:
        t=templates[r['grating']]; wins=[(5000,5450),(5600,5850),(6000,6400)]
        if r['grating']==300:wins.append((4500,4800))
        a={f'{lo}_{hi}':shift_fit(r['w'],r['y'],t['w'],t['y'],lo,hi) for lo,hi in wins}
        for v in a.values():v['conditional_barycentric_relative_kms']=v['shift_kms']+r['bcv']-t['bcv']
        q={k:v for k,v in r.items() if k not in ['w','y','flux']};q.update(template=t['file'],windows=a)
        if r['grating']==300:q['telluric_6860_6910']=shift_fit(r['w'],r['y'],t['w'],t['y'],6860,6910)
        q['median_shift']=float(np.median([v['shift_kms'] for v in a.values()]))
        q['window_range']=float(np.ptp([v['shift_kms'] for v in a.values()]))
        q['conditional_barycentric_median']=q['median_shift']+r['bcv']-t['bcv']
        records.append(q)
    inject=[]
    for g,t in templates.items():
        fake=np.interp(t['w']/(1+10/C),t['w'],t['y'])
        for lo,hi in [(5000,5450),(5600,5850),(6000,6400)]:inject.append(dict(grating=g,window=[lo,hi],injected=10,**shift_fit(t['w'],fake,t['w'],t['y'],lo,hi)))
    sensitivity=[]
    for r in rows:
        if r['grating']!=300:continue
        for width in [40,70,100]:
            t=next(x for x in rows if x['file'].startswith('0146.'))
            shifts=[shift_fit(r['w'],normalized(r['w'],r['flux'],width),t['w'],normalized(t['w'],t['flux'],width),lo,hi)['shift_kms'] for lo,hi in [(5000,5450),(5600,5850),(6000,6400),(4500,4800)]]
            sensitivity.append(dict(file=r['file'],reference=t['file'],continuum_A=width,shifts=shifts,median=float(np.median(shifts)),spread=float(np.ptp(shifts))))
    blurred=[]
    for r in rows:
        if r['grating']!=300:continue
        t=templates[300]
        yy=gaussian_filter1d(r['y'],7/2.35482/np.median(np.diff(r['w'])))
        tt=gaussian_filter1d(t['y'],7/2.35482/np.median(np.diff(t['w'])))
        shifts=[shift_fit(r['w'],yy,t['w'],tt,lo,hi)['shift_kms'] for lo,hi in [(5000,5450),(5600,5850),(6000,6400),(4500,4800)]]
        blurred.append(dict(file=r['file'],shifts=shifts,median=float(np.median(shifts)),spread=float(np.ptp(shifts))))
    result=dict(qualification='Apparent relative spectral shifts only. Conditional barycentric column assumes native observed-frame wavelength; absolute zero point, slit-centering and LSF are uncalibrated. No orbit fit or astrophysical RV error claimed.',records=records,injections=inject,sensitivity_direct_repeat_template=sensitivity,sensitivity_extra_FWHM7A=blurred)
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in records:print(r['file'],r['date'],r['grating'],'median',round(r['median_shift'],2),'bary',round(r['conditional_barycentric_median'],2),'range',round(r['window_range'],2),'sky',round(r['sky']['center_A'],3),'tell',round(r.get('telluric_6860_6910',{}).get('shift_kms',0),2))
if __name__=='__main__':main()

