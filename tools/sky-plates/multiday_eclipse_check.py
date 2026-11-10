"""Existing extracted TESS multi-day eclipse sensitivity, not completeness/FAP."""
import json
from pathlib import Path
import numpy as np
from analyze_tess_variability import load_products, bin_series, bin_values

P = 35.9204566
REF = 2500.
DUR = [1., 2., 4., 6.]
PHASES = np.arange(288) * P / 288

def trapezoid(t, center, duration, period=P):
    d = abs((np.asarray(t)-REF-center+period/2) % period-period/2)
    return np.clip((duration/2-d)/(.2*duration), 0, 1)

def projection(t, group, degree):
    cols=[]
    for g in np.unique(group):
        m=group==g; x=(t-np.mean(t[m]))/10
        cols.extend(m*x**k for k in range(degree+1))
    x=np.column_stack(cols)
    u,s,_=np.linalg.svd(x,full_matrices=False)
    return u[:,s>max(x.shape)*np.finfo(float).eps*s[0]]

def residualize(a,q):
    return a-q@(q.T@a)

def run(products, degree, period=P):
    phasegrid=np.arange(288)*period/288
    ts=[];ys=[];groups=[];banks=[]
    for i,p in enumerate(products):
        b=bin_series(p['native_t'],p['native_y'],1/24)
        ts.extend(b['t']);ys.extend(b['y']);groups.extend([i]*len(b['t']))
        bank=np.column_stack([bin_values(trapezoid(p['native_t'],ph,d,period),b['inverse'],len(b['t'])) for d in DUR for ph in phasegrid])
        banks.append(bank)
    t=np.array(ts);y=np.array(ys);g=np.array(groups);bank=np.vstack(banks)
    q=projection(t,g,degree); sr=residualize(-bank,q); yr=residualize(y,q)
    norms=np.sum(sr*sr,axis=0); valid=norms>1e-8
    coeff=np.divide(sr.T@yr,norms,out=np.zeros(len(norms)),where=valid)
    score=coeff*np.sqrt(norms); best=int(np.argmax(score))
    durations=np.repeat(DUR,len(phasegrid)); phases=np.tile(phasegrid,len(DUR))
    observed=dict(phase_days=float(phases[best]),duration_days=float(durations[best]),depth=float(coeff[best]),score=float(score[best]),baseline_rms=float(np.sqrt(np.mean(yr**2))))
    inject_idx=np.array([j*len(PHASES)+k for j in range(len(DUR)) for k in range(0,len(PHASES),12)])
    cross=sr.T@sr[:,inject_idx]
    invnorm=np.divide(1,norms,out=np.zeros(len(norms)),where=valid)
    depth_grid=[.001,.003,.01,.03,.1];out=[]
    for depth in depth_grid:
        injected_coeff=coeff[:,None]+depth*cross*invnorm[:,None]
        injscore=injected_coeff*np.sqrt(norms[:,None])
        win=np.argmax(injscore,axis=0)
        for j,d in enumerate(DUR):
            sl=slice(j*24,(j+1)*24);w=win[sl];true=inject_idx[sl]
            distance=abs((phases[w]-phases[true]+period/2)%period-period/2)
            recovered=injected_coeff[w,np.arange(j*24,(j+1)*24)]
            success=(distance<d/2)&(durations[w]>=d/2)&(durations[w]<=d*2)&(recovered>=depth*.5)&(recovered<=depth*2)&(injscore[w,np.arange(j*24,(j+1)*24)]>score[best])
            out.append(dict(duration_days=d,injected_depth=depth,recovered_phases=int(sum(success)),tested_phases=24,successful_phase_days=phases[true][success].tolist(),depth_at_true_template_minmax=[float(np.min(injected_coeff[true,np.arange(j*24,(j+1)*24)])),float(np.max(injected_coeff[true,np.arange(j*24,(j+1)*24)]))]))
    coverage=[]
    for j,d in enumerate(DUR):
        sub=bank[:,j*len(PHASES):(j+1)*len(PHASES)]
        coverage.append(dict(duration_days=d,occupied_hour_bins_in_eclipse_minmax=[int(np.min(np.sum(sub>0,axis=0))),int(np.max(np.sum(sub>0,axis=0)))],projected_template_energy_minmax=[float(min(norms[j*288:(j+1)*288])),float(max(norms[j*288:(j+1)*288]))]))
    return dict(products=[p['key'] for p in products],period_days=period,degree=degree,n_hour_bins=len(t),observed_best=observed,injections=out,coverage=coverage,observation_phases=[dict(duration_days=d,phase_days=float(phases[j*288+k]),depth=float(coeff[j*288+k]),score=float(score[j*288+k])) for j,d in enumerate(DUR) for k in range(288)])

def main():
    root=Path('results/sky-plates-j0754-multiday-eclipses-20261109');root.mkdir(exist_ok=True)
    products,inventory=load_products(Path('results/sky-plates-j0754-checks-20261108'))
    outputs=[]
    for provider in ['QLP','TESS-SPOC']:
        for corrected in [False,True]:
            a=[p for p in products if p['provider']==provider and p['primary']==corrected]
            for degree in [1,2]:
                r=run(a,degree);r.update(provider=provider,corrected=corrected);outputs.append(r)
                print(provider,corrected,degree,r['observed_best'],flush=True)
    (root/'results.json').write_text(json.dumps(outputs,indent=2,allow_nan=False))
    (root/'native_inventory.json').write_text(json.dumps(inventory,indent=2))
    a=[p for p in products if p['provider']=='QLP' and not p['primary']]
    sensitivity=[run(a,degree,period) for period in [P-.030613065,P+.030613065] for degree in [1,2]]
    (root/'period_sensitivity.json').write_text(json.dumps(sensitivity,indent=2,allow_nan=False))
    # Post-hoc stability diagnostic at the selected pooled QLP SAP trough.
    selected=outputs[0]['observed_best']; per_sector=[]
    for p in products:
        b=bin_series(p['native_t'],p['native_y'],1/24)
        s=-bin_values(trapezoid(p['native_t'],selected['phase_days'],selected['duration_days']),b['inverse'],len(b['t']))
        for degree in [1,2]:
            q=projection(b['t'],np.zeros(len(b['t'])),degree)
            sr=residualize(s,q); yr=residualize(b['y'],q); energy=sr@sr
            per_sector.append(dict(key=p['key'],degree=degree,energy=float(energy),depth=float(sr@yr/energy) if energy>1e-8 else None))
    (root/'selected_trough_per_sector.json').write_text(json.dumps(per_sector,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for ax,base in zip(axes.ravel(),range(0,8,2)):
        for r in outputs[base:base+2]:
            a=[v for v in r['observation_phases'] if v['duration_days']==6]
            ax.plot([v['phase_days'] for v in a],[100*v['depth'] for v in a],label=f"sector degree {r['degree']}")
        ax.axhline(0,c='gray',lw=.5)
        ax.set_title(outputs[base]['provider']+(' corrected' if outputs[base]['corrected'] else ' SAP'))
        ax.set_xlabel('Mid-eclipse phase (days from BTJD 2500)')
        ax.set_ylabel('Fitted 6-day depth (%)');ax.legend()
    fig.suptitle('Descriptive template depths: reduction and baseline sensitivity, not eclipse detections')
    fig.savefig(root/'six_day_profiles.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
