"""Bounded existing-data V/g and optical/IR trajectory diagnostics, no network."""
import csv,json,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from asassn_replication import parse_row

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-chromaticity-20261109'
P=35.9204566

def design(t,camera,band):
    t=np.asarray(t);camera=np.asarray(camera);band=np.asarray(band)
    columns=[]
    for c in sorted(set(camera)):
        m=camera==c;columns.extend([m.astype(float),m*(t-np.mean(t[m]))/180])
    for b in ['V','g']:
        m=band==b;phase=2*np.pi*(t-2457000)/P
        columns.extend([m*np.cos(phase),m*np.sin(phase)])
    return np.column_stack(columns)

def fit_optical(path):
    allrows=[parse_row(r) for r in csv.DictReader(path.open())]
    accepted=[r for r in allrows if r['accepted']]
    selected=[r for r in accepted if '2017-10-07'<=r['day']<'2018-05-01']
    groups=defaultdict(list)
    for r in selected:groups[(r['camera'],r['day'])].append(r)
    daily=[dict(t=np.mean([r['jd'] for r in a]),y=np.median([r['mag'] for r in a]),camera=c,band=a[0]['band'],day=d) for (c,d),a in groups.items()]
    counts=Counter(r['camera'] for r in daily)
    daily=[r for r in daily if counts[r['camera']]>=10]
    t,y=np.array([[r['t'],r['y']] for r in daily]).T
    X=design(t,[r['camera'] for r in daily],[r['band'] for r in daily]);op=np.linalg.pinv(X)
    coef=op@y;res=y-X@coef
    out={'range':['2017-10-07','2018-05-01'],'daily_counts':dict(Counter((r['camera']+'_'+r['band']) for r in daily)),
         'coefficients_V_g':coef[-4:].tolist(),'amplitudes_V_g':np.linalg.norm(coef[-4:].reshape(2,2),axis=1).tolist(),
         'g_minus_V_coefficients':(coef[-2:]-coef[-4:-2]).tolist(),'block_sensitivity':{}}
    rng=np.random.default_rng(20261109)
    for days in [7,14]:
        blocks=np.floor((t-t.min())/days).astype(int);draws=[]
        for _ in range(1000):
            signs=rng.choice([-1,1],blocks.max()+1)
            draws.append((op@(X@coef+res*signs[blocks]))[-4:])
        draws=np.array(draws);diff=draws[:,2:]-draws[:,:2]
        out['block_sensitivity'][str(days)]={'amplitude_V_g_q025_q975':np.quantile(np.linalg.norm(draws.reshape(-1,2,2),axis=2),[.025,.975],axis=0).tolist(),
            'g_minus_V_coeff_q025_q975':np.quantile(diff,[.025,.975],axis=0).tolist(),
            'qualification':'Shared-time block-sign residual sensitivity, not calibrated confidence intervals.'}
    injected=np.zeros(X.shape[1]);injected[-4:]=[.03,-.01,.05,-.02]
    out['exact_injection_coefficient_max_error']=float(np.max(np.abs(op@(X@injected)-injected)))
    out['rank']=int(np.linalg.matrix_rank(X));out['columns']=X.shape[1]
    v=[r for r in accepted if r['band']=='V'];g=[r for r in accepted if r['band']=='g']
    pairs=[(r,min(g,key=lambda z:abs(z['jd']-r['jd']))) for r in v]
    close=[(a,b) for a,b in pairs if abs(a['jd']-b['jd'])<=.5]
    out['overlap_inventory']={'V_dates':[min(r['day'] for r in v),max(r['day'] for r in v)],'g_dates':[min(r['day'] for r in g),max(r['day'] for r in g)],'nearest_pairs_within_half_day':len(close),'distinct_V_days':len(set(a['day'] for a,b in close)),'pairs_within_0p05day':sum(abs(a['jd']-b['jd'])<=.05 for a,b in close)}
    return out

def predictions(g,g0,g1,ir0,ir1):
    x=10**(-.4*(np.asarray(g)-g0));x1=10**(-.4*(g1-g0));y1=10**(-.4*(ir1-ir0))
    affine=1+(y1-1)*(x-1)/(x1-1)
    return ir0-2.5*np.log10(affine),ir0+(ir1-ir0)*(np.asarray(g)-g0)/(g1-g0)

def trajectory(path):
    d=json.loads(path.read_text());out=[]
    for window,values in d['optical_matches'].items():
        for camera in sorted(set(r['camera'] for r in values['matches'])):
            rows=[r for r in values['matches'] if r['camera']==camera]
            a=[r for r in rows if r['year']==2018];b=[r for r in rows if r['year']==2022]
            if not a or not b:continue
            g0,g1=[float(np.median([r['mag'] for r in z])) for z in [a,b]]
            i0,i1=[float(np.median([r['infrared_mag'] for r in z])) for z in [a,b]]
            held=[r for r in rows if r['year'] not in [2018,2022]]
            pred,linear=predictions([r['mag'] for r in held],g0,g1,i0,i1)
            residual=np.array([r['infrared_mag'] for r in held])-pred
            out.append({'window_days':int(window),'camera':camera,'training_years':[2018,2022],'training_g':[g0,g1],'training_W1':[i0,i1],
                'n_evaluation':len(held),'flux_affine_rms_mag':float(np.sqrt(np.mean(residual**2))),
                'mag_linear_rms_mag':float(np.sqrt(np.mean((np.array([r['infrared_mag'] for r in held])-linear)**2))),
                'max_prediction_difference_mag':float(np.max(np.abs(pred-linear))),
                'rows':[{'date':r['date'][:10],'g':r['mag'],'W1':r['infrared_mag'],'affine_prediction':float(p),'residual':float(e)} for r,p,e in zip(held,pred,residual)]})
    return out

def main():
    ap=ROOT/'results/sky-plates-j0754-independent-20261108/asassn_v1_lightcurve.csv'
    dp=ROOT/'results/sky-plates-j0754-dust-20261108/dust_analysis.json'
    zp=ROOT/'results/sky-plates-j0754-checks-20261108/ztf_667401913089980544.csv'
    z=list(csv.DictReader(zp.open()))
    result={'optical':fit_optical(ap),'trajectory':trajectory(dp),'ztf_filter_counts':dict(Counter(r['filtercode'] for r in z)),
        'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ap,dp,zp]}}
    OUT.mkdir(exist_ok=True);(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
