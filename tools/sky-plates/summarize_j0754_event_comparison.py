"""Readable four-event report and conditional, profile-specific rate summary."""
import hashlib,json
from pathlib import Path
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/sky-plates-j0754-flare-comparison-20261111'

def main():
    r=json.loads((OUT/'comparison.json').read_text());lines=[
        'J0754: COMMON MEASUREMENTS OF FOUR TESS EXCURSIONS',
        '11 November 2026 | Existing calibrated pixels | Exploratory follow-up','',
        'RESULT',
        'All four selected excursions retain positive, compact excesses. November22',
        'and December23 have the clearest fast-rise/slow-decay profiles. November16',
        'is broader and nearly symmetric; November28 has a moderate rise/decay',
        'asymmetry. The largest event remains December23. If these are stellar flares',
        'associated with the J0754 system, their conditional radiated energies are',
        'in the superflare regime. The observations do not resolve which component',
        'of a possible close binary emits. No new star or companion is identified.','',
        'COMMON MODEL MEASUREMENTS',
        'Date         Peak     FWHM     Rise tau  Decay tau   Model ED   Observed ED',
        '             percent  hours    hours     hours       seconds    seconds']
    for e in r['events']:
        f=e['photometry']['central'];lines.append(f"{e['date_utc'][:10]}   {100*f['amplitude_fractional']:5.2f}   {f['fwhm_h']:5.2f}     {f['rise_tau_h']:5.2f}     {f['decay_tau_h']:5.2f}       {f['equivalent_duration_s']:6.1f}      {f['observed_signed_ED_s']:6.1f}")
    lines+=['','Model ED integrates the fitted two-sided exponential to its full tails.',
        'Observed ED is a signed integral over -12..+18 hours relative to the selected',
        'peak, excluding gaps longer than2.5native cadences. It is neither clipped',
        'positive nor filled through missing intervals. November22 lacks much of the',
        'earlier context; December23 contains a50-minute retained-sample gap.',
        'Reported FWHM is the intrinsic model width, not a separately measured',
        'physical flare duration. The model is a shape description, not a proof of',
        'magnetic reconnection. Native frame integration is included in the fit.','',
        'UNCERTAINTY: BASELINE/APERTURE SENSITIVITY IS MATERIAL',
        'Date         Conditional bootstrap ED16..84%    All tested model ED range']
    for e in r['events']:
        b=e['photometry']['bootstrap']['p16_p50_p84'];lo,hi=e['energy']['model_ED_baseline_aperture_range_s']
        lines.append(f"{e['date_utc'][:10]}              {b[0][2]:.0f}..{b[2][2]:.0f}s                        {lo:.0f}..{hi:.0f}s")
    lines+=['','The bootstrap resamples one-hour residual blocks at a fixed model/baseline.',
        'It is not the total uncertainty. The variant envelope uses three aperture',
        'radii(1,1.5,2pixels), linear/quadratic baselines and two context windows.',
        'All48models converge without boundary solutions and retain positive areas.',
        'The broad November16 and November28 areas can roughly double under these',
        'choices; November22 and December23 are more stable. These are sensitivity',
        'ranges, not confidence intervals. No selection-adjusted false-alarm',
        'probability is assigned.','',
        'CONDITIONAL BOLOMETRIC RADIATED ENERGY',
        'Date         Face-on surface patch     Isotropic equivalent',
        '             erg                       erg']
    for e in r['events']:
        n=e['energy'];lines.append(f"{e['date_utc'][:10]}       {n['face_on_patch_erg']:.3e}                 {n['isotropic_equivalent_erg']:.3e}")
    lines+=['','Nominal assumptions shared by every event: Lstar54.640694Lsun,',
        'Tstar5205.7764K,Tflare9000K and the retained SVO TESS energy response.',
        'E_iso=ED*Lstar*(Tflare/Tstar)^4*band(Tstar)/band(Tflare).',
        'The face-on optically thick surface-patch convention is E_iso/4; this is',
        'a geometry convention, not four independently measured energies. Varying',
        'L30..60Lsun,Tstar4500..5250K,Tflare8000..12000K and baseline/aperture gives',
        'the event-specific envelopes in comparison.json. Unknown dilution,',
        'epoch-dependent stellar luminosity, flare spectrum and geometry remain',
        'systematic assumptions. This measures radiation, not total magnetic or',
        'kinetic energy. The smallest tested face-on energy is about9e36erg.',
        'The older December23 ED~601s/patch energy~5.3e37erg used a finite window',
        'and different local baseline. Its nominal value remains inside the new',
        'systematic envelope. The new full-model ED~809s/patch energy~7.1e37erg',
        'must not be described as an additional event or an independent energy test.','',
        'POSITION AND POINTING CHECK',
        'Date         Bias-corrected template offset   Matched-injection95% radius']
    for e in r['events']:
        l=e['localization'];lines.append(f"{e['date_utc'][:10]}              {l['calibrated_offset_pixels']:.4f}pixels                     {l['injection_radius95_pixels']:.4f}pixels")
    lines+=['','Offsets are relative to the observed target-image template. Each event has',
        '80blank-time injected source-shaped controls with exactly the same retained',
        'before/event/after counts and time pattern. All finite native frames must',
        'match within0.3cadence; flagged, missing and duplicated matches are rejected.',
        'This is empirical-template calibration, not an official TESS PRF posterior.',
        'The December23 offset is marginally outside its noise-only radius, by',
        'about0.0006pixel; template/color/astrometric systematics are not included.',
        'The November28 localization is appreciably less precise. None of these',
        'positions resolves a close stellar companion.','',
        'HOURLY DISPLACEMENT IS NOT FULLY SETTLED',
        'Date         Speed(pixels/hour)    Conditional stationary-reference p']
    for e in r['events']:
        m=e['localization']['hourly_motion'];lines.append(f"{e['date_utc'][:10]}             {m['speed_pixels_per_h']:.4f}                       {m['empirical_stationary_reference_p']:.4f}")
    lines+=['','November22 and December23 show apparent hourly shifts unusual relative to',
        'the matched stationary controls. These discrete empirical reference values',
        'are not calibrated significance levels or proof of a moving celestial',
        'object. The hourly positions also have non-linear scatter. Image-gradient',
        'terms representing small pointing/scene changes alter event-mean source',
        'amplitudes by less than0.3% relative; all selected hourly changes are below',
        '0.9%. Thus the tested pointing terms do not remove the flux excess, while',
        'stationary-source astrometry is not established by this check. A formal',
        'PRF/attitude analysis or moving-object ephemeris exclusion remains open.',
        'The new result strengthens photometric comparison, not a claim that every',
        'possible contaminant has been excluded.','',
        'DETECTION SENSITIVITY AND EXPOSURE',
        'The direct-pixel QUALITY=0 series contains114.51nominal accepted-cadence',
        'days across six disjoint intervals. This differs from104.18days in the',
        'earlier union of extracted light curves because their quality selections',
        'are different. Neither number is detector live-time. QLP/SPOC alternatives',
        'are reductions of the same exposures, not extra observing nights.',
        'Four measured profiles were injected at48uniform peak times in each of',
        'six native direct-aperture series:1152trials before median removal.','',
        'Injected profile    Recovery(all times)   Recovery(>=80% area sampled)   Effective days']
    exposures={}
    for e in r['events']:
        rows=[x for x in r['sensitivity']['rows'] if x['event_date']==e['date_utc'][:10]]
        nr=sum(x['new_recoveries'] for x in rows);nt=sum(x['trials'] for x in rows);ns=sum(x['recovered_well_sampled'] for x in rows);nw=sum(x['well_sampled_trials'] for x in rows);te=sum(x['effective_search_days'] for x in rows);exposures[e['date_utc'][:10]]=te
        lines.append(f"{e['date_utc'][:10]}                 {nr}/{nt}                  {ns}/{nw}                        {te:.2f}")
    te=exposures['2021-12-23'];rate=1/te;interval=[chi2.ppf(.025,2)/(2*te),chi2.ppf(.975,4)/(2*te)]
    rate_summary=dict(scope='Illustration assuming one accepted December23-like event and a constant Poisson occurrence rate for copies of that exact measured profile.',count=1,effective_search_days=te,rate_per_day=rate,conditional_poisson95_per_day=interval,
        qualification='Finite-grid recovery exposure, candidate attribution and fixed-profile assumptions are not included in these Poisson intervals. Not a fitted flare-frequency distribution or measured magnetic-cycle dependence.')
    (OUT/'conditional_rate.json').write_text(json.dumps(rate_summary,indent=2)+'\n',encoding='utf8',newline='\n')
    lines+=['','These direct-pixel trials are separate from the1296earlier extracted-',
        'light-curve injections. Recovery is specific to these shapes and this',
        'screen. The broad/baseline-sensitive event models can be easier to recover',
        'than other real flares with the same energy. No upstream FFI calibration',
        'transfer or complete flare-frequency distribution is inferred.','',
        f'For repeated copies of the December23 profile, one accepted event over',
        f'{te:.2f}effective search days gives an illustrative rate {rate:.5f}/day.',
        f'The conditional Poisson95% interval is {interval[0]:.5f}..{interval[1]:.5f}/day.',
        'This assumes a constant rate and correct attribution; it is too sparse to',
        'measure a magnetic-cycle correlation. Four candidate dates divided by',
        'coverage alone would ignore differences in profile and detection sensitivity.','',
        'REPRODUCTION',
        'python -B tools/sky-plates/compare_j0754_tess_events.py --sensitivity',
        'python -B tools/sky-plates/summarize_j0754_event_comparison.py',
        'python -B -m unittest discover -s tools/sky-plates -p test_*.py','',
        'All native pixels remain unchanged. comparison.json records exact inputs,',
        'parameters, variants, bootstrap draws, injections and source/protocol hashes.',
        'Recovery trials were reused unchanged while correcting the position fitter;',
        'the reuse chain is retained. --reuse-sensitivity is for deliberate unchanged-',
        'profile reruns; a fresh reproduction should use --sensitivity.',
        'Figures: four_event_measurements.png; four_event_localization.png.','',
        'PRIMARY METHOD REFERENCES',
        'https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowdingv2.html',
        'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=TESS/TESS.Red',
        'https://doi.org/10.1051/0004-6361/202449405',
        'https://doi.org/10.1051/0004-6361/202040098','']
    (OUT/'report.txt').write_text('\n'.join(lines),encoding='utf8',newline='\n')
    verification=dict(comparison_source_matches=hashlib.sha256((ROOT/'tools/sky-plates/compare_j0754_tess_events.py').read_bytes()).hexdigest()==r['source_sha256'],
        exact_control_counts=all(all(x['matched_counts']==e['localization']['observed_before_event_after_counts'] for x in e['localization']['controls']) for e in r['events']),
        original_input_hashes_match=all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in r['input_sha256'].items()),
        bandpass_sha256=hashlib.sha256((ROOT/'results/sky-plates-j0754-flare-energy-20261108/TESS_Red.xml').read_bytes()).hexdigest(),
        summary_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        event_count=len(r['events']),fit_variant_count=sum(len(e['photometry']['variants']) for e in r['events']))
    if not all(verification[k] for k in ['comparison_source_matches','exact_control_counts','original_input_hashes_match']):raise ValueError('Verification mismatch')
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(verification=verification,rate=rate_summary),indent=2))

if __name__=='__main__':main()
