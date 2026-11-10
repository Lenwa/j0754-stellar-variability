# J0754: measured radial-velocity information in the existing FAST spectra

9 November 2026. This is a new relative-line-shift diagnostic using the ten already saved spectra, not an orbit fit or a request for new data.

## What the data can determine

The spectra do contain measurable velocity-sensitive structure. Numerical injections of a 10 km/s shift are recovered to numerical precision in both grating groups. Low instrumental resolution alone is therefore not a sufficient reason to declare that the spectra contain no velocity information.

However, actual same-night repeatability and cross-window consistency are inadequate for a calibrated test of the candidate 10.72 km/s orbital semi-amplitude. This is an empirical result of examining the spectra, rather than an assumed limitation based only on their resolution.

## Native inputs and method

The ten original FAST FITS remain unchanged and their SHA256 hashes are in results.json. The 300-line/mm and 600-line/mm groups are fitted separately. The linear native wavelength solution and extracted spectrum are used. Continuum removal divides by a broad Gaussian-smoothed spectrum (70 Angstrom standard deviation), followed by fitted multiplicative line amplitude and additive offset in each comparison window. An unweighted least-squares shift searches -150 to +150 km/s, with a coarse grid followed by local refinement.

Photospheric windows are 5000–5450, 5600–5850 and 6000–6400 Angstrom; the 300-line group also uses 4500–4800 Angstrom. These are diagnostic blended-line windows, not calibrated stellar templates. The high-SNR 0135 spectrum is the initial 300-line reference and 0089 the 600-line reference. No H-alpha velocity, archive galaxy classification or automated scalar pipeline velocity is adopted as astrophysical truth.

The headers contain BCV and XCSAO output, including template-relative velocities and corrected reported velocities. BCV differences are retained as an explicitly conditional observed-frame correction; the wavelength processing frame itself is not independently certified by this audit. The decisive same-night result is effectively independent of that uncertainty because the BCV change is only about 0.018 km/s.

## Same-night repeatability: a direct calibration test

Spectra 0146, 0147 and 0148 were taken on 12 November 2009, spanning about 15.3 minutes in exposure midpoint. Their apparent median shifts relative to the older reference are -55.4, -67.6 and -79.4 km/s, a range of about 24.0 km/s. Disjoint-window spreads are approximately 42.3, 43.4 and 29.2 km/s.

Directly comparing the repeat spectra to 0146 rather than the older template still gives apparent median shifts of -23.3 km/s (0147) and -20.3 km/s (0148). The latter's disjoint-window range is about 7.8 km/s; the former's is 48.1 km/s. Varying continuum smoothing from 40 to 100 Angstrom changes these medians by only about 0.2 km/s. Thus simple continuum choice or the older reference is not enough to resolve the inconsistency.

For a circular 35.9205-day, K=10.7245 km/s orbit, the maximum change over this 15.3-minute interval is approximately 0.020 km/s. The observed apparent shift differences cannot plausibly be the velocity change of that proposed orbit. They expose noise/extraction/slit illumination/line-shape or wavelength-systematic limitations in this use of the data. They do not diagnose which effect dominates.

## Across nights and grating groups

Most non-reference 300-line spectra have disjoint-window ranges of approximately 29–92 km/s. Adding an extra Gaussian FWHM of 7 Angstrom alters some medians substantially and leaves large inconsistencies (including a 116 km/s window range in one 2012 spectrum). This is a sensitivity check, not a claim to have measured or matched the true line-spread functions.

One promising pair is the 600-line January 2010 pair: 0159 relative to 0089 has an apparent shift of 5.38 km/s, or 7.08 km/s after the conditional differential BCV correction, with only 3.22 km/s spread across the three windows. But the other 600-line exposure has a 60.9 km/s window range. A coherent pair is insufficient to certify the zero point or establish an orbit; the 7.08 km/s value is not a calibrated stellar RV measurement.

## Sky and telluric checks

The extracted background contains a measurable 5577-Angstrom sky feature. Gaussian fits have peak amplitude/residual-RMS ratios of about 25–107, with no fitted-center boundary failures. The same-night sky centers move only about 0.12 Angstrom, or 6.6 km/s, while the stellar diagnostics move much more. Sky emission fills the slit differently from a stellar point source, so a clean sky correction cannot remove unknown stellar centering errors.

For the low-resolution spectra, a 6860–6910 Angstrom telluric-region comparison gives same-night shifts approximately -33.3, -51.4 and -43.7 km/s against the old reference. Subtracting these from the photospheric medians still leaves a repeat range of roughly 20 km/s. This narrow region is a blended diagnostic, not a verified isolated atmospheric fiducial; its wavelength-dependent line profile and stellar blends preclude claiming that it supplies a calibrated correction.

## Result and limits

These existing spectra can establish that the naive relative-RV method fails repeatability at the scale needed for this particular orbital question. They can also identify the January 2010 higher-resolution pair as the cleanest retained relative-shift comparison for future reprocessing of the same data. They do not independently confirm or reject Gaia's orbit, measure a trustworthy stellar RV amplitude, or distinguish orbital motion from spots.

No phase extrapolation was attempted. Gaia's correlated near-circular periastron parameters cannot be propagated independently, and its period uncertainty also matters over the earlier FAST years. No formal RV error is quoted: synthetic same-template shift recovery tests numerical interpolation, not real extraction, slit centering, template mismatch or calibration errors.

## Saved state

Protocol, full shifts, window results, same-night/continuum/smoothing sensitivity, sky diagnostics, metadata and original input hashes are in results.json and protocol.txt in this directory. Maintained source is tools/sky-plates/fast_relative_rv.py, with test_fast_relative_rv.py and fast_relative_rv_README.md. Three new tests cover shift recovery, amplitude/offset invariance and sky-center fitting. No external retrieval, contact, original-data modification or publication occurred.

Verification: all 37 photospheric and 7 telluric-region winning template amplitudes are positive, so no anticorrelation minimum is being reported as a shift. Parent source review found no blocking issue. The current 135-test component suite passed, including three new FAST tests. Scoped source commit 30de3bcb7645cfe640609f56430df01df7eae951; own files clean, unrelated edits preserved, no push. Source hashes and verification are saved in verification.json.
