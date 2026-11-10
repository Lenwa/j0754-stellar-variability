"""Existing-data exploratory oscillation information; no seismic classification."""
from pathlib import Path
import json, hashlib
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-oscillation-check-20261109'
INPUT=ROOT/'results/sky-plates-j0754-checks-20261108'
FREQ=np.arange(10.,200.001,.05)

def residual(t,y):
    x=np.column_stack([np.ones(len(t)),t-np.mean(t)])
    return y-x@np.linalg.lstsq(x,y,rcond=None)[0]

def comb(t,rms=1e-4):
    rng=np.random.default_rng(8413)
    freqs=np.arange(24.,97.,6.)
    amp=np.exp(-.5*((freqs-60)/15)**2)
    y=np.sum(amp[:,None]*np.cos(2*np.pi*freqs[:,None]*1e-6*(t-t.min())[None,:]*86400+rng.uniform(0,2*np.pi,len(freqs))[:,None]),axis=0)
    return rms*y/np.std(y)

def spectrum(t,y):
    # PSD-like scaling: periodogram power times twice effective cadence.
    # Fractional flux squared -> ppm squared is 1e12; per Hz -> per microHz
    # is 1e-6. Net factor 1e6 below. Only relative band shape/ratios are
    # interpreted for irregular cadence; no calibrated absolute PSD inference.
    return LombScargle(t-t.min(),y,normalization='psd').power(FREQ*1e-6*86400,method='fast')*2*np.ptp(t)*86400/len(t)*1e6

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    records=[]; curves=[]; hashes={}
    for path in sorted(INPUT.glob('*lc.fits')):
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as h:
            d=h[1].data; sector=h[0].header.get('SECTOR',h[1].header.get('SECTOR'))
            if sector is None:
                sector=int(path.name.split('_s')[1][:4])
            source='QLP' if 'qlp' in path.name else 'SPOC'
            for col in ['SAP_FLUX','KSPSAP_FLUX','DET_FLUX','PDCSAP_FLUX']:
                if col not in d.names:continue
                rawt=np.asarray(d['TIME'],float); rawf=np.asarray(d[col],float)
                mask=(d['QUALITY']==0)&np.isfinite(rawt)&np.isfinite(rawf)&(rawf>0)
                base=mask.copy();mask &= ~((rawt>=2571.5)&(rawt<=2573.0))
                t=rawt[mask];y=rawf[mask]/np.median(rawf[mask])-1
                order=np.argsort(t);t=t[order];y=y[order]
                yr=residual(t,y);power=spectrum(t,yr)
                smooth=gaussian_filter1d(power,5/.05)
                selection=(FREQ>=25)&(FREQ<=110)
                peak=float(FREQ[selection][np.argmax(smooth[selection])])
                bands={f'{lo}-{lo+20}':float(np.mean(power[(FREQ>=lo)&(FREQ<lo+20)])) for lo in range(20,181,20)}
                inj=[]
                for rms in [1e-4,3e-4,1e-3]:
                    signal=comb(t,rms);processed=residual(t,signal)
                    sigpower=spectrum(t,processed)
                    s=(FREQ>=30)&(FREQ<=100)
                    inj.append(dict(rms_ppm=rms*1e6,linear_detrend_rms_retention=float(np.std(processed)/np.std(signal)),injected_only_to_observed_30_100_power=float(np.sum(sigpower[s])/np.sum(power[s]))))
                rec=dict(source=source,sector=int(sector),column=col,n=len(t),span_days=float(np.ptp(t)),median_cadence_minutes=float(np.median(np.diff(t))*1440),residual_rms_ppm=float(np.std(yr)*1e6),difference_scatter_ppm=float(np.median(abs(np.diff(yr)-np.median(np.diff(yr))))*1.4826/np.sqrt(2)*1e6),smoothed_peak_25_110_uHz=peak,bands=bands,injections=inj)
                if np.any(base!=mask):
                    tb=rawt[base];yb=rawf[base]/np.median(rawf[base])-1
                    pb=spectrum(tb,residual(tb,yb));s=(FREQ>=30)&(FREQ<=100)
                    rec['flare_retained_to_masked_band_power']=float(np.sum(pb[s])/np.sum(power[s]))
                records.append(rec);curves.append((source,sector,col,power,smooth))
    (OUT/'results.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    (OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    np.savez_compressed(OUT/'power_spectra.npz',frequency_uHz=FREQ,**{f'{a}_{b}_{c}':p for a,b,c,p,s in curves})
    sectors=sorted(set(r['sector'] for r in records));fig,axes=plt.subplots(3,2,figsize=(13,11),layout='constrained')
    for ax,sector in zip(axes.flat,sectors):
        for source,sec,col,p,s in curves:
            if sec==sector:ax.semilogy(FREQ,s,label=f'{source} {col}',lw=1)
        ax.set_title(f'Sector {sector}');ax.set_xlim(10,200);ax.set_xlabel('Frequency (microHz)');ax.set_ylabel('Smoothed PSD-like power');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Existing TESS oscillation-shape diagnostic; separate reductions, no detection threshold')
    fig.savefig(OUT/'power_shapes.png',dpi=140);plt.close(fig)
    print(json.dumps([{k:r[k] for k in ['source','sector','column','smoothed_peak_25_110_uHz','residual_rms_ppm']} for r in records],indent=2))

if __name__=='__main__':main()
