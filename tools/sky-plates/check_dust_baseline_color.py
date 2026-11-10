"""Conditional absolute-color check on the high-baseline dust rescue."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse
from stellar_models import band_flux,read_bandpass
from compare_stellar_models import dust_grid,extinction_deficits

def fnu(wave,flux,band):
    # A constant F_nu spectrum has F_lambda = c F_nu / lambda^2.
    return band_flux(wave,flux,band)/band_flux(wave,2.99792458e18/wave**2,band)

def main(root,out):
    models=root/'results/sky-plates-j0754-models-20261108/models'
    g=read_bandpass(models/'filter_g.xml','SLOAN/SDSS.g');w=read_bandpass(models/'filter_W1.xml','WISE/WISE.W1')
    v=parse(models/'filter_W1.xml');params={p.name:p.value for p in v.iter_fields_and_params() if hasattr(p,'value')}
    assert params['MagSys']=='Vega' and params['ZeroPointUnit']=='Jy'
    zp=float(params['ZeroPoint']);r=json.loads((out/'results.json').read_text())
    chosen=[x for x in r['variable_column_trajectories'] if x['camera']=='bq' and x['window']==30 and x['tau0']==3]
    rows=[]
    for t in [4500,5000,5500]:
        with fits.open(models/f'ckm05_{t}.fits') as h:wave=np.array(h[1].data['WAVELENGTH'],float);s=np.array(h[1].data['g25'],float)
        intrinsic=-2.5*np.log10(fnu(wave,s,g)/fnu(wave,s,w))+2.5*np.log10(3631/zp)
        for fit in chosen:
            grain=dust_grid(root/f"results/sky-plates-j0754-dust-20261108/draine_{fit['law']}_extinction.txt")
            d=extinction_deficits(wave,s,[g,w],grain,np.array([3]))[0];tr=1-fit['cover']*d
            color=float(intrinsic-2.5*np.log10(tr[0]/tr[1]))
            rows.append(dict(T=t,law=fit['law'],cover=fit['cover'],intrinsic_gAB_W1Vega=float(intrinsic),predicted_color=color,observed_color=12.621-9.50675,residual=color-(12.621-9.50675)))
    result=dict(W1_zero_Jy=zp,g_AB_zero_Jy=3631,rows=rows,
        caveat='Nominal5000K temporal fits tested against baseline color;4500/5500K sensitivity holds geometry fixed, not refitted. ASASSNg uses SDSSg proxy and nominalAB convention; using savedSVOgVega zero4023.573 instead adds0.1115mag to modelcolor. No foreground reddening, extra emission or secondstar. Source photometric offsets/response differences not fully calibrated; not formalrejection.',
        calibration_source='https://www.astronomy.ohio-state.edu/asassn/public/examples.shtml',
        input_hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [out/'results.json',models/'filter_g.xml',models/'filter_W1.xml',*[models/f'ckm05_{t}.fits' for t in [4500,5000,5500]]]})
    (out/'baseline_color.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.root,a.out)
