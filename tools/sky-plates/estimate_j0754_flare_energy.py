"""Conditional blackbody radiated energy from the preserved optical event."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io.votable import parse_single_table

OUT = Path('results/sky-plates-j0754-flare-energy-20261108')
INPUT = Path('results/sky-plates-j0754-physical-tests-20261108/tess/sector46_pixel_event_timeseries.csv')


def main():
    OUT.mkdir(exist_ok=True)
    data = np.genfromtxt(INPUT, delimiter=',', names=True, dtype=None, encoding='utf8')
    t = (data['btjd']-2572.0782528500076)*24
    y = data['relative_core_flux']
    pre = (t >= -12) & (t <= -4)
    post = (t >= 14) & (t <= 17)
    # Two simple local baselines plus the existing pixel reduction's baseline.
    slope = (np.median(y[post])-np.median(y[pre]))/(np.median(t[post])-np.median(t[pre]))
    baselines = {'preserved_reference': np.zeros_like(t),
        'local_pre_constant': np.full_like(t,np.median(y[pre])),
        'local_pre_post_linear': np.median(y[pre])+slope*(t-np.median(t[pre]))}
    tab = parse_single_table(str(OUT/'TESS_Red.xml')).array
    wave = np.asarray(tab['Wavelength'],float)*1e-10
    response = np.asarray(tab['Transmission'],float)
    def band(T):
        # SVO declares DetectorType=0 (energy response); do not add lambda.
        planck = wave**-5/np.expm1(.014387768775039337/(wave*T))
        return np.trapezoid(response*planck,wave)
    def correction(Tstar,Tflare):
        return (Tflare/Tstar)**4*band(Tstar)/band(Tflare)
    records = []
    for name,b in baselines.items():
        for start,stop in [(-3,12),(-3,14),(-3,16)]:
            use = (t>=start)&(t<=stop)
            # No positive clipping; retain noise to avoid upward integration bias.
            ed = float(np.trapezoid((y[use]-b[use])/(1+b[use]),t[use]*3600))
            maxgap = float(np.max(np.diff(t[use])))
            assert ed>0 and maxgap<1
            records.append(dict(baseline=name,window_hours=[start,stop],equivalent_duration_seconds=ed,
                                maximum_gap_hours=maxgap,gap_handling='Trapezoid linearly bridges retained-sample gaps, including one ~50-minute gap.'))
    central = next(r for r in records if r['baseline']=='local_pre_post_linear' and r['window_hours']==[-3,14])
    L = 54.640694*3.828e33
    Tstar,Tflare = 5205.7764,9000
    energy = central['equivalent_duration_seconds']*L*correction(Tstar,Tflare)
    grid = [r['equivalent_duration_seconds']*lum*3.828e33*correction(ts,tf)
        for r in records for lum in (30,54.640694,60)
        for ts in (4500,5205.7764,5250) for tf in (8000,9000,10000,12000)]
    result = dict(status='Approximate isotropic-equivalent bolometric radiated energy; conditional on event belonging to J0754 and assumed blackbody spectra.',
        isotropic_equivalent_energy_erg=energy,isotropic_equivalent_energy_joules=energy*1e-7,
        projected_surface_patch_energy_erg=energy/4,
        geometry='Isotropic equivalent uses 4pi d^2. Standard face-on optically thick flare-patch convention uses projected area pi Rstar^2 and A sigma Tflare^4, giving one quarter of this. Unknown surface orientation changes inferred emitted energy.',
        equivalent_duration_seconds=central['equivalent_duration_seconds'],
        central_assumptions=dict(stellar_luminosity_Lsun=54.640694,stellar_temperature_K=Tstar,flare_temperature_K=Tflare,
            stellar_provenance='Preserved Gaia DR3 astrophysical_parameters FLAME luminosity and GSP-Phot temperature; systematic uncertainty exceeds formal catalogue intervals.'),
        correction_factor=correction(Tstar,Tflare),
        isotropic_sensitivity_envelope_erg=[min(grid),max(grid)],
        surface_patch_sensitivity_envelope_erg=[min(grid)/4,max(grid)/4],sensitivity_grid='L=30..60 Lsun; Tstar=4500..5250 K; Tflare=8000..12000 K; nine baseline/window combinations. Illustrative systematic envelope, not confidence interval.',
        integrations=records,formula='E_iso = ED * Lstar_bol * (Tflare/Tstar)^4 * integral(R Bstar)/integral(R Bflare); E_face_on_patch=E_iso/4',
        source_input=str(INPUT),input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        response_url='https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=TESS/TESS.Red',
        method_reference='https://doi.org/10.1051/0004-6361/202142088',
        geometry_reference='https://doi.org/10.1051/0004-6361/202243869',response_detector_type=0,
        limits='No new image reduction, extinction correction, aperture/PRF throughput correction, neighbour exclusion, or statistical flare-template fit. Assumes relative aperture variation represents whole-star fraction; different flare/stellar aperture throughput and unknown spectral evolution add systematic error. This estimates radiation, not kinetic or total magnetic energy.')
    (OUT/'energy_estimate.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
