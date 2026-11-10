"""Cadence/annual-alias sensitivity; no period significance or orbital confirmation."""
import csv, hashlib, json
from pathlib import Path
import numpy as np
from astropy.time import Time
from orbit_photometry import nuisance, daily

ROOT=Path('results/sky-plates-j0754-alias-check-20261109')
P=35.9204566
PERIODS={'Gaia':P,'year_plus':1/(1/P+1/365.25),'year_minus':1/(1/P-1/365.25),'half':P/2,'double':2*P,'lunar':29.53059}

def orth(a):
    u,s,_=np.linalg.svd(a,full_matrices=False)
    return u[:,:np.sum(s>max(a.shape)*np.finfo(float).eps*s[0])]

def prepare(t,e,g):
    return orth(nuisance(t,g,'group_linear')/e[:,None])

def residual(y,e,q):
    z=y/e
    return z-q@(q.T@z)

def basis(t,e,q,p):
    angle=2*np.pi*(t-2457000)/p
    a=np.column_stack([np.cos(angle),np.sin(angle)])/e[:,None]
    return orth(a-q@(q.T@a))

def power(z,b):
    return float(np.sum((b.T@z)**2)/(z@z))

def run(rows):
    t=np.array([float(r['t']) for r in rows]); y=np.array([float(r['y']) for r in rows])
    e=np.maximum(.02,[float(r['err']) for r in rows]); g=np.array([r['group'] for r in rows])
    dates=Time(t,format='jd').isot
    season=np.array([int(d[:4])-(int(d[5:7])<7) for d in dates])
    records=[]; leak=[]; rng=np.random.default_rng(20261109)
    for year in sorted(set(season)):
        m=season==year
        if m.sum()<30 or np.ptp(t[m])<120: continue
        tt,yy,ee,gg=t[m],y[m],e[m],g[m]
        q=prepare(tt,ee,gg); z=residual(yy,ee,q)
        bs=[basis(tt,ee,q,p) for p in PERIODS.values()]
        scores={name:power(z,b) for name,b in zip(PERIODS,bs)}
        win={name:float(abs(np.mean(np.exp(2j*np.pi*(tt-tt.min())/p)))**2) for name,p in {'annual':365.25,'lunar':29.53059,'Gaia':P}.items()}
        matrix=np.zeros((16,len(PERIODS),len(PERIODS)))
        for j,p in enumerate(PERIODS.values()):
            for k in range(16):
                phase=rng.uniform(0,2*np.pi)
                injected=np.cos(2*np.pi*(tt-2457000)/p+phase)
                zz=residual(injected,ee,q)
                matrix[k,j,:]=[power(zz,b) for b in bs]
        leak.append(matrix)
        records.append({'season':str(year),'n':int(m.sum()),'span_days':float(np.ptp(tt)),'power':scores,'window_power':win})
    avg=np.mean(leak,axis=0)
    margins=np.array([r['power']['Gaia']-max(r['power']['year_plus'],r['power']['year_minus']) for r in records])
    return {'gaia_beats_both_annual_count':int(np.sum(margins>0)),'leave_one_out_min_mean_margin':float(min(np.mean(np.delete(margins,i)) for i in range(len(margins)))),'seasons':records,'mean_power':{name:float(np.mean([r['power'][name] for r in records])) for name in PERIODS},'injection_mean_crosstalk':np.mean(avg,axis=0).tolist(),'injection_phase_min':np.min(avg,axis=0).tolist(),'injection_phase_max':np.max(avg,axis=0).tolist(),'qualification':'Rows injected periods, columns fitted periods; equal-season mean explained fraction of nuisance-projected noiseless sinusoid. Phase ranges are 16 chosen simulations, not confidence intervals.'}

def main():
    a=Path('results/sky-plates-j0754-models-20261108/short_period/asassn_g_daily.csv')
    z=Path('results/sky-plates-j0754-checks-20261108/modern_differential.json')
    with a.open() as f: rows=list(csv.DictReader(f))
    raw=json.loads(z.read_text())['measurements']
    zrows=daily([dict(t=r['mjd']+2400000.5,y=r['corrected_relative_mag'],err=.02,group='ZTF') for r in raw])
    result={'periods_days':PERIODS,'ASASSN':run(rows),'ZTF':run(zrows),'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a,z]}}
    ROOT.mkdir(exist_ok=True)
    (ROOT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]['mean_power'] for k in ['ASASSN','ZTF']},indent=2))
if __name__=='__main__': main()

