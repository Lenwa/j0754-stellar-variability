"""Passband integration and bounded spectral-change models; no causal priors."""
from dataclasses import dataclass
import numpy as np
from astropy.io.votable import parse


@dataclass
class Bandpass:
    name: str
    wave: np.ndarray
    transmission: np.ndarray
    detector: int  # SVO: 0 energy, 1 photon; spectra use F_lambda.


def read_bandpass(path,expected_id):
    v=parse(path)
    params={x.name:x.value for x in v.iter_fields_and_params() if hasattr(x,'value')}
    if params.get('filterID')!=expected_id:
        raise ValueError(f'Wrong response identity: {path}: {params.get("filterID")}')
    if str(params.get('WavelengthUnit')).lower()!='angstrom':
        raise ValueError('Expected Angstrom response wavelength')
    detector=int(params['DetectorType'])
    if detector not in [0,1]:raise ValueError('Unknown detector weighting')
    t=v.get_first_table().to_table()
    return Bandpass(expected_id,np.asarray(t['Wavelength'],float),np.asarray(t['Transmission'],float),detector)


def band_flux(wave,flux,band):
    wave=np.asarray(wave,float);flux=np.asarray(flux,float)
    if np.any(np.diff(wave)<=0) or len(wave)!=len(flux):raise ValueError('Invalid spectrum grid')
    positive=band.transmission>0
    if not np.any(positive) or wave[0]>min(band.wave[positive]) or wave[-1]<max(band.wave[positive]):
        raise ValueError('Spectrum does not cover response')
    x=np.unique(np.r_[band.wave,wave[(wave>=band.wave.min())&(wave<=band.wave.max())]])
    response=np.interp(x,band.wave,band.transmission,left=0,right=0)
    weight=response*(x if band.detector==1 else 1)
    result=np.trapezoid(np.interp(x,wave,flux)*weight,x)/np.trapezoid(weight,x)
    if not np.isfinite(result) or result<=0:raise ValueError('Nonpositive band flux')
    return float(result)


def magnitudes(faint,bright):
    faint=np.asarray(faint,float);bright=np.asarray(bright,float)
    if np.any(faint<=0) or np.any(bright<=0):raise ValueError('Positive flux required')
    return -2.5*np.log10(faint/bright)


def interpolate(temperatures,fluxes,temperature):
    if temperature<temperatures[0] or temperature>temperatures[-1]:raise ValueError('Temperature outside grid')
    return np.array([np.interp(temperature,temperatures,v) for v in np.asarray(fluxes).T])


def fit_temperature_scale(temperatures,fluxes,bright_temperature,observed):
    bright=interpolate(temperatures,fluxes,bright_temperature)
    desired=bright[0]/bright[1]*10**(-.4*(observed[0]-observed[1]))
    residual=fluxes[:,0]-desired*fluxes[:,1]
    roots=[]
    for i in range(len(temperatures)-1):
        if temperatures[i]>bright_temperature:break
        if residual[i]==0:roots.append(float(temperatures[i]))
        elif residual[i]*residual[i+1]<0:
            t=temperatures[i]-residual[i]*(temperatures[i+1]-temperatures[i])/(residual[i+1]-residual[i])
            if t<=bright_temperature:roots.append(float(t))
    if not roots:return None
    faint_temperature=max(roots)
    unscaled=interpolate(temperatures,fluxes,faint_temperature)
    delta=magnitudes(unscaled,bright)
    gray=float(np.mean(np.asarray(observed)[:2]-delta[:2]))
    return dict(temperature_faint=faint_temperature,gray_magnitude=gray,
                area_ratio=float(10**(-.4*gray)),prediction=delta+gray,solutions=len(roots))


def fit_partial_cover(tau,deficit,observed):
    loss=1-10**(-.4*np.asarray(observed))
    best=None
    for t,d in zip(tau,deficit):
        if d[0]<=0:continue
        fraction=loss[0]/d[0]
        if not 0<=fraction<=1:continue
        prediction=-2.5*np.log10(1-fraction*d)
        score=float(np.sum((prediction[:2]-observed[:2])**2))
        if best is None or score<best['score']:
            best=dict(tau=float(t),covering_fraction=float(fraction),prediction=prediction,score=score)
    return best


def fit_two_component(hot,cool,observed):
    q=np.asarray(cool)/np.asarray(hot)
    loss=1-10**(-.4*np.asarray(observed))
    denominator=loss[0]*q[0]-loss[1]*q[1]
    if abs(denominator)<1e-14:return None
    area=(loss[1]-loss[0])/denominator
    cover=loss[0]*(1+area*q[0])
    if area<0 or not 0<=cover<=1:return None
    prediction=-2.5*np.log10(1-cover/(1+area*q))
    return dict(companion_area_ratio=float(area),covered_fraction=float(cover),prediction=prediction)


def fit_spots(hot,cool_grid,temperatures,observed,start=0):
    loss=1-10**(-.4*np.asarray(observed));best=None
    for t,cool in zip(temperatures,cool_grid):
        deficit=1-np.asarray(cool)/np.asarray(hot)
        if deficit[0]<=0 or np.any(1-start*deficit<=0):continue
        change=loss[0]*(1-start*deficit[0])/deficit[0]
        end=start+change
        if not start<=end<=1:continue
        prediction=-2.5*np.log10((1-end*deficit)/(1-start*deficit))
        score=float(np.sum((prediction[:2]-observed[:2])**2))
        if best is None or score<best['score']:
            best=dict(spot_temperature=float(t),spot_fraction_bright=start,spot_fraction_faint=float(end),prediction=prediction,score=score)
    return best


def covering_predictions(deficits,observed,tolerance=.02,limit=1.):
    """All sampled opacity/contrast cases with physical covering intervals.

    At each supplied spectral deficit retain both bounds and the midpoint of
    the fraction interval allowed by g and W1, including non-optimal solutions.
    No check-band value enters these bounds. This is grid/tolerance sensitivity,
    not a statistical confidence set or a completeness claim between grid nodes.
    """
    observed=np.asarray(observed)
    low_loss=1-10**(-.4*(observed[:2]-tolerance))
    high_loss=1-10**(-.4*(observed[:2]+tolerance))
    result=[]
    for i,d in enumerate(deficits):
        if np.any(d[:2]<=0):continue
        lo=max(0.,float(np.max(low_loss/d[:2])))
        hi=min(limit,float(np.min(high_loss/d[:2])))
        if lo>hi:continue
        for fraction in [lo,(lo+hi)/2,hi]:
            faint=1-fraction*d
            if np.any(faint<=0):continue
            result.append(dict(index=i,fraction=fraction,prediction=-2.5*np.log10(faint)))
    return result
