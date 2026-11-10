"""Conditional physical equations tied to saved J0754 observations.

This module calculates hypotheses, not orbit validation or confidence limits.
The CLI reuses declared g/W1 fits and extends them with Cousins R predictions.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy import constants as const
from astropy.io import fits
from scipy.optimize import brentq

from stellar_models import band_flux, interpolate, read_bandpass
from compare_stellar_models import dust_grid, extinction_deficits

G=const.G.value
MSUN=const.M_sun.value
RSUN=const.R_sun.value
AU=const.au.value
DAY=86400.


def pressure_scale_height(temperature,logg_cgs,mu=1.3):
    """Ideal-gas, local hydrostatic pressure scale height, metres."""
    if temperature<=0 or mu<=0:raise ValueError('Positive temperature and molecular weight required')
    return float(const.k_B.value*temperature/(mu*const.u.value*(10**logg_cgs/100)))


def radius_from_logg(mass_solar,logg_cgs):
    if mass_solar<=0:raise ValueError('Positive mass required')
    return float(np.sqrt(G*mass_solar*MSUN/(10**logg_cgs/100))/RSUN)


def logg_from_radius(mass_solar,radius_solar):
    if mass_solar<=0 or radius_solar<=0:raise ValueError('Positive mass and radius required')
    return float(np.log10(G*mass_solar*MSUN/(radius_solar*RSUN)**2*100))


def ellipsoidal_coefficient(limb_u,gravity_g):
    """Leading fractional FLUX SEMI-amplitude coefficient; not 15/32."""
    if not 0<=limb_u<=1 or gravity_g<0:raise ValueError('Unsupported darkening coefficients')
    return float(.15*(15+limb_u)*(1+gravity_g)/(3-limb_u))


def mass_function(period_days,K_kms,eccentricity):
    if period_days<=0 or K_kms<0 or not 0<=eccentricity<1:raise ValueError('Invalid orbit inputs')
    return float(period_days*DAY*(K_kms*1000)**3*(1-eccentricity**2)**1.5/(2*np.pi*G*MSUN))


def minimum_companion_mass(primary_mass,f_mass):
    if primary_mass<=0 or f_mass<=0:raise ValueError('Positive mass and mass function required')
    return float(brentq(lambda m:m**3/(primary_mass+m)**2-f_mass,1e-8,max(100.,100*f_mass)))


def roche_fraction(primary_to_companion_mass):
    """Eggleton circular/synchronous volume radius of the PRIMARY / separation."""
    if primary_to_companion_mass<=0:raise ValueError('Positive mass ratio required')
    q=primary_to_companion_mass**(1/3)
    return float(.49*q*q/(.6*q*q+np.log1p(q)))


def inclination_from_K(observed_K,edge_K):
    if observed_K<0 or edge_K<=0:raise ValueError('Invalid RV amplitudes')
    if observed_K>edge_K*(1+1e-12):return None
    return float(np.degrees(np.arcsin(min(1.,observed_K/edge_K))))


def planck_contrast(wavelength_um,cool_temperature,hot_temperature):
    w=np.asarray(wavelength_um,float)
    if np.any(w<=0) or min(cool_temperature,hot_temperature)<=0:raise ValueError('Positive wavelengths and temperatures required')
    return np.expm1(14387.768775/(w*hot_temperature))/np.expm1(14387.768775/(w*cool_temperature))


def spot_prediction(contrasts,reference_fade):
    contrasts=np.asarray(contrasts,float)
    deficit=1-contrasts
    fraction=(1-10**(-.4*reference_fade))/deficit[0]
    if not 0<=fraction<=1:raise ValueError('Spot coverage cannot reproduce anchor')
    return dict(spot_fraction=float(fraction),fade_mag=(-2.5*np.log10(1-fraction*deficit)).tolist())


def partial_powerlaw(fade_a,fade_b,wavelength_a,wavelength_b,covering_fraction):
    loss=1-10**(-.4*np.array([fade_a,fade_b]))
    if not 0<covering_fraction<=1 or np.any(loss>=covering_fraction) or np.any(loss<=0):
        raise ValueError('Finite positive depths require covering greater than each flux loss')
    tau=-np.log1p(-loss/covering_fraction)
    alpha=-np.log(tau[1]/tau[0])/np.log(wavelength_b/wavelength_a)
    return dict(covering_fraction=float(covering_fraction),tau_reference=float(tau[0]),
                reference_wavelength_um=float(wavelength_a),alpha=float(alpha))


def partial_prediction(wavelengths_um,model):
    tau=model['tau_reference']*(np.asarray(wavelengths_um)/model['reference_wavelength_um'])**(-model['alpha'])
    return -2.5*np.log10(1-model['covering_fraction']*(-np.expm1(-tau)))


def photometric_radius(Vmag,distance_pc,Av,BCv,temperature):
    """Conditional V/BC/extinction/temperature calculation, not an observed radius."""
    if distance_pc<=0 or temperature<=0:raise ValueError('Positive distance and temperature required')
    mv=Vmag-5*np.log10(distance_pc/10)-Av
    luminosity=10**((4.74-(mv+BCv))/2.5)
    return dict(absolute_V=float(mv),absolute_bolometric=float(mv+BCv),
                luminosity_solar=float(luminosity),radius_solar=float(np.sqrt(luminosity)*(5772/temperature)**2))


def binary_case(m1,m2,radius_solar,period_days,K_kms,eccentricity,limb_u=.8,gravity_g=.32):
    separation=(G*(m1+m2)*MSUN*(period_days*DAY)**2/(4*np.pi**2))**(1/3)
    edgeK=2*np.pi*separation/(period_days*DAY)*m2/(m1+m2)/np.sqrt(1-eccentricity**2)/1000
    inclination=inclination_from_K(K_kms,edgeK)
    roche=separation/RSUN*roche_fraction(m1/m2)
    alpha=ellipsoidal_coefficient(limb_u,gravity_g)
    edge_ell=alpha*m2/m1*(radius_solar*RSUN/separation)**3
    ell=None if inclination is None else edge_ell*np.sin(np.radians(inclination))**2
    contact_ratio=(radius_solar+.5)*RSUN/separation
    return dict(M1_solar=m1,M2_solar=m2,radius_solar=radius_solar,separation_AU=float(separation/AU),
                separation_Rsun=float(separation/RSUN),primary_roche_Rsun=float(roche),
                roche_filling=float(radius_solar/roche),K_edge_kms=float(edgeK),inclination_deg=inclination,
                ellipsoidal_alpha=alpha,ellipsoidal_fractional_flux_semiamplitude_edge=float(edge_ell),
                ellipsoidal_fractional_flux_semiamplitude_at_candidate_i=None if ell is None else float(ell),
                ellipsoidal_approximation_caution=bool(radius_solar/roche>.8),
                central_circular_first_fourth_contact_days=float(period_days/np.pi*np.arcsin(contact_ratio)) if contact_ratio<1 else None,
                eclipse_min_inclination_deg=float(np.degrees(np.arccos(contact_ratio))) if contact_ratio<1 else None,
                candidate_i_eclipse_possible=bool(inclination is not None and contact_ratio<1 and np.cos(np.radians(inclination))<=contact_ratio))


def extend_R_predictions(model_report,model_dir,dust_dir,Rband):
    """Evaluate the SAME prior fitted hypotheses in a new response; no W2 fitting."""
    cache={};grains={};records=[]
    temperatures=np.arange(3500.,6501.,250.)
    for source_index,old in enumerate(model_report['records']):
        if old['response']!='g' or old['season']!='March':continue
        key=(old['metallicity'],old['logg'])
        if key not in cache:
            spectra=[]
            for t in temperatures:
                with fits.open(model_dir/f'{key[0]}_{int(t)}.fits') as h:
                    wave=np.asarray(h[1].data['WAVELENGTH'],float)
                    spectra.append(np.asarray(h[1].data[key[1]],float))
            grid=np.array([band_flux(wave,s,Rband) for s in spectra])
            cache[key]=(wave,np.array(spectra),grid)
        wave,spectra,grid=cache[key];hotT=old['temperature_bright'];hot=interpolate(temperatures,grid[:,None],hotT)[0]
        p=old['parameters'];family=old['family']
        if family in ['temperature_plus_scale','fixed_area_cooling']:
            cool=interpolate(temperatures,grid[:,None],p['temperature_faint'])[0]
            delta=-2.5*np.log10(cool/hot)+p.get('gray_magnitude',0.)
        elif family=='cool_spots':
            cool=interpolate(temperatures,grid[:,None],p['spot_temperature'])[0];contrast=cool/hot
            delta=-2.5*np.log10((1-p['spot_fraction_faint']*(1-contrast))/(1-p['spot_fraction_bright']*(1-contrast)))
        elif family=='two_component_occultation':
            cool=interpolate(temperatures,grid[:,None],old['companion_temperature'])[0]
            delta=-2.5*np.log10(1-p['covered_fraction']/(1+p['companion_area_ratio']*cool/hot))
        elif family in ['partial_dust_covering','uniform_dust_screen']:
            law=old['dust_law']
            if law not in grains:grains[law]=dust_grid(dust_dir/f'draine_{law}_extinction.txt')
            spec=interpolate(temperatures,spectra,hotT)
            deficit=extinction_deficits(wave,spec,[Rband],grains[law],np.array([p['tau']]))[0,0]
            delta=-2.5*np.log10(1-p.get('covering_fraction',1.)*deficit)
        elif family=='gray_single_component':delta=(.312+.1175)/2
        else:raise ValueError('Unknown existing model family '+family)
        bands=dict(old['bands'],R=float(delta))
        records.append(dict(family=family,source_record_index=source_index,metallicity=key[0],logg=key[1],temperature_bright=hotT,
                            **{k:old[k] for k in ['dust_law','companion_temperature'] if k in old},
                            parameters=p,bands=bands,training_residual_max=old['max_training_residual'],
                            W2_residual=old['W2_residual'],g_minus_R=bands['g']-bands['R'],R_minus_I=bands['R']-bands['I']))
    summary={}
    for family in sorted({r['family'] for r in records}):
        compatible=[r for r in records if r['family']==family and r['training_residual_max']<=.02]
        summary[family]=dict(cases=len([r for r in records if r['family']==family]),g_W1_compatible=len(compatible),
            predictions={band:[min(r['bands'][band] for r in compatible),max(r['bands'][band] for r in compatible)] for band in ['g','R','r','I','W1','W2']} if compatible else {},
            R_minus_I=[min(r['R_minus_I'] for r in compatible),max(r['R_minus_I'] for r in compatible)] if compatible else None)
    return records,summary


def main(models,dust,model_report,literature,out):
    out.mkdir(parents=True,exist_ok=True)
    audit=json.loads(literature.read_text());gaia=audit['nss_quality'];context=audit['Gaia_source_context']
    observed={'g':.312,'W1':.1175,'W2':.08875};waves={'g':.48,'V':.55,'R':.65,'I':.8,'W1':3.4,'W2':4.6}
    contrast=planck_contrast(list(waves.values()),3800.,4500.);spot=spot_prediction(contrast,observed['g'])
    spot.update(temperature_phot_K=4500.,temperature_spot_K=3800.,contrast=dict(zip(waves,contrast.tolist())),
                fade_mag=dict(zip(waves,spot['fade_mag'])),qualification='Monochromatic illustration anchored on actual g, not an observed V fade or fitted atmosphere.')
    partial=[]
    for fraction in [1.,.5,.3,.26]:
        m=partial_powerlaw(observed['g'],observed['W1'],waves['g'],waves['W1'],fraction)
        m['fade_mag']=dict(zip(waves,partial_prediction(list(waves.values()),m).tolist()));partial.append(m)
    radius_cases=[]
    for mass in [.8,1.2,2.6]:
        for label,radius,temperature in [('poster_logg1.7_conditional_mass',radius_from_logg(mass,1.7),4500.),('Gaia_GSPphot_cautioned',audit['atmosphere_catalogue']['radius_gspphot'],audit['atmosphere_catalogue']['teff_gspphot'])]:
            logg=logg_from_radius(mass,radius);height=pressure_scale_height(temperature,logg)
            radius_cases.append(dict(label=label,mass_solar=mass,radius_solar=radius,temperature_K=temperature,logg_cgs=logg,
                scaleheight_m=height,scaleheight_Rsun=height/RSUN,H_over_R=height/(radius*RSUN),illustrative_5_to_10_H_over_R=[5*height/(radius*RSUN),10*height/(radius*RSUN)]))
    fmass=mass_function(gaia['period'],gaia['semi_amplitude_primary'],gaia['eccentricity'])
    binary=[];minimum=[]
    for mass in [.8,1.2,2.6]:
        minmass=minimum_companion_mass(mass,fmass);minimum.append(dict(M1_solar=mass,minimum_M2_solar=minmass))
        for m2 in [minmass,.5,1.]:
            for label,radius in [('poster_logg1.7',radius_from_logg(mass,1.7)),('Gaia_GSPphot_cautioned',audit['atmosphere_catalogue']['radius_gspphot'])]:
                binary.append(dict(radius_assumption=label,**binary_case(mass,m2,radius,gaia['period'],gaia['semi_amplitude_primary'],gaia['eccentricity'])))
    distance=1000/context['parallax'];photo=[]
    for temperature in [4500.,5200.]:
        for Av in [0.,.15,.5]:
            for BCv in [-.7,-.5,-.3]:
                V=audit['Tang2010_constraints']['V_ASAS_median']
                photo.append(dict(Vmag=V,distance_pc=distance,Av_assumed=Av,BCv_assumed=BCv,temperature_assumed=temperature,**photometric_radius(V,distance,Av,BCv,temperature)))
    Rband=read_bandpass(out/'filter_CousinsR.xml','Generic/Cousins.R')
    records,summary=extend_R_predictions(json.loads(model_report.read_text()),models,dust,Rband)
    report=dict(observed_fades_mag=observed,observed_faint_to_bright_flux={k:10**(-.4*v) for k,v in observed.items()},
        planck_fixed_4500_3800=spot,partial_powerlaw_illustrations=partial,integrated_atmosphere_R_extension=summary,
        atmosphere_records=records,scaleheight_cases=radius_cases,photometric_radius_sensitivity=photo,
        Argon_assumed_V12p35_reproduction=photometric_radius(12.35,1668.,.15,-.5,4500.),
        candidate_Gaia_NSS=gaia,candidate_mass_function_solar=fmass,minimum_companion_masses=minimum,binary_cases=binary,
        ellipsoidal_prefactor_audit=dict(correct='.15*(15+u)*(1+g)/(3-u)',supplied='15*(15+u)*(1+g)/(32*(3-u))',supplied_over_correct=3.125),
        limitations=['No actual V=.35 measurement was used; g approximated by .48um only in Planck/power-law illustrations, integrated grid retains SDSSg response proxy.',
            'R means Generic/Cousins.R; Sloanr is a separate existing passband. R forecasts reuse prior fits without W2 selection; ranges are per-nuisance best fits, not confidence sets or full tolerance envelopes.',
            'Gaia NSS remains a candidate, and historical8km/s RV precision is not K<=8. Mass and radius assumptions are not observed posterior distributions.',
            'Ellipsoidal leading term is fractional flux semi-amplitude, assumes small tidal distortion and undiluted primary; near-Roche cases violate this approximation and cannot be quoted as exact photometry.',
            'Scale height uses local hydrostatic ideal gas and assumed mu;5-10 heights is illustrative. Optical depth does not convert directly into pressure/radius without opacity/atmosphere modeling; logg discrepancy does not diagnose NLTE.',
            'V_ASAS12.24 is a historical median, not simultaneous with Gaia or modern fading; inverse parallax, BC,Av,T assumptions make radius sensitivity conditional.',
            'Spot, partial-cover and temperature/scale explanations remain degenerate. A free area factor is not a measured radius change. No unique extinction exponent implies a unique grain size.',
            'No SPHEREx refit or dust mass/detection is performed here.'],
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [model_report,literature,out/'filter_CousinsR.xml']})
    (out/'physical_model_results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(spot=spot,partial=partial,mass_function=fmass,minimum_companions=minimum,R_forecasts=summary),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--models',type=Path,required=True);parser.add_argument('--dust',type=Path,required=True)
    parser.add_argument('--model-report',type=Path,required=True);parser.add_argument('--literature',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();main(args.models,args.dust,args.model_report,args.literature,args.out)
