# J0754 physical diagnostics on preserved observations

For a simple optical flare shape comparison, run
`compare_optical_flare_shapes.py` from repository root with the same Python runtime.
It reads the preserved sector-46 pixel time series and writes to
`results/sky-plates-j0754-optical-overlay-20261108`. The solar curve is an
approximate manual digitization of Hao et al. (2017), Fig. 3c; the stellar
curve uses Davenport et al. (2014), equations 1 and 4. Peak/FWHM alignment
compares shapes only; it is neither an energy estimate nor a flare confirmation.

`estimate_j0754_flare_energy.py` integrates that same preserved time series
with three baseline choices and three integration windows. It reads the SVO
`TESS_Red.xml` response in `results/sky-plates-j0754-flare-energy-20261108`
and writes a conditional blackbody energy receipt there. It distinguishes
isotropic-equivalent radiation from the face-on surface-patch convention;
its sensitivity envelope is not a statistical confidence interval.

`check_fast_bluewing.py` checks native FAST Halpha profiles by grating, registers
nearby metal lines, and compares independent-night profiles against other nights.
The exploratory protocol and outputs live in
`results/sky-plates-j0754-fast-bluewing-20261108`. Residuals are descriptive:
this does not supply calibrated outflow velocities or detection limits.

This extends the existing sky-plate investigation with the tests Kyle requested
after supplying Argon's mathematical audit. The audit's tables and limits are
external claims; none are substituted for local measurements. Maintained source
is this component. Results, protocols, native inputs, models and Python runtime
remain outside source Git in the repository root.

Use the existing `python` with this component's
`requirements.txt`. No additional dependency installation is required. These
commands run from repository root; all input snapshots must already exist.

```powershell
$skyPython = 'python'
$skySource = 'tools/sky-plates'
$skyPhysical = 'results/sky-plates-j0754-physical-tests-20261108'
& $skyPython -B "$skySource/fast_line_tests.py" --input 'results/sky-plates-j0754-spherex-20261108/literature/fast' --out "$skyPhysical/fast"
& $skyPython -B "$skySource/analyze_tess_variability.py" --out "$skyPhysical/tess"
& $skyPython -B "$skySource/analyze_tess_variability.py" --out "$skyPhysical/tess" --pixels-only "$skyPhysical/tess/pixels/tess-s0046-1-4_118.691174_16.694621_7x7_astrocut.fits"
& $skyPython -B "$skySource/analyze_spherex_co.py" --previous 'results/sky-plates-j0754-spherex-20261108' --models 'results/sky-plates-j0754-models-20261108/models' --dust 'results/sky-plates-j0754-dust-20261108' --out "$skyPhysical/co"
& $skyPython -B "$skySource/physical_constraints.py" --models 'results/sky-plates-j0754-models-20261108/models' --dust 'results/sky-plates-j0754-dust-20261108' --model-report 'results/sky-plates-j0754-models-20261108/model_comparison.json' --literature 'results/sky-plates-j0754-dust-20261108/literature/j0754_literature_audit.json' --out "$skyPhysical/models"
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

FAST files contain one-dimensional extracted spectra, raw/background/sigma
bands and a linear wavelength calibration. TESS files contain one-dimensional
QLP or TESS-SPOC light curves, with quality/background/centroid diagnostics.
They are telescope measurements, not a fresh reduction of two-dimensional
detector pixels. The ten FAST spectra represent eight nights; paired TESS
reductions share photons. Neither duplication creates independent epochs.

A targeted follow-up also retrieved a public 7x7 calibrated FFI pixel cube for
sector46. The optional `--pixels-only` command uses that preserved cube and a
pre-inspection protocol amendment to compare the fixed positive-event window
with guarded adjacent references. Difference-image moments, background-plane
and motion-derivative sensitivities support only conditional localization. They
do not provide PRF-calibrated astrometry, a confirmed flare, or raw detector counts.

`fast_line_tests.py` matches approximate instrumental resolution, propagates
the full smoothing/resampling covariance and fitted-continuum gradient into
equivalent-width errors, and reports independent-night contrasts. Native sigma
independence and approximate line-spread functions remain assumptions. Numerical
injections distinguish broad additive background, local emission and continuum
placement errors. Sideband, resolution, wavelength and background variants are
retained. Broad blended indices do not establish precise temperatures, resolved
line-depth ratios or radial velocities.

`analyze_tess_variability.py` retains sector/reduction identity, screens native
quality without clipping real large positive/negative excursions, and preserves
sampling gaps. Joint signal/nuisance fits use the existing candidate-period basis.
Event/injection diagnostics describe response of extracted light curves, not
upstream pipeline completeness, global false-alarm probabilities, or a universal
transit bound. Equivalent durations are observed-aperture quantities, not
bolometric flare energies. No long-term DC baseline is stitched across sectors.
The observed-residual bootstrap includes candidate excursions; none above its
screen cannot establish flare absence. One-day high-pass and 1/3/6-hour box
injections do not calibrate multi-day giant-binary eclipse sensitivity.

`analyze_spherex_co.py` extends the previous static near-IR atmosphere sensitivity
to logg1–3 and evaluates a 4.55–4.9 micron/reference4.0–4.2 micron residual index.
Current and native wavelength labels/reductions remain separate. Approximate
rectangular responses, coarse LTE atmospheres, foreground extinction and model
degeneracy prevent identifying CO emission or measuring dust mass. A 5% near-IR
shape tolerance is a grid sensitivity set, not a confidence region. This test
does not reconstruct the missing 2022 spectrum.

`physical_constraints.py` evaluates scale heights, conditional radii, Roche
geometry, candidate mass function and ellipsoidal flux semi-amplitudes. It uses
actual g/W1/W2 changes and the candidate catalogue K, not the attachment's
unverified V fade or RV upper bound. The ellipsoidal approximation requires a
detached, modestly distorted system; Roche-filling cases cannot validate it.
Pressure scale height is not a measured atmospheric extension. Model R/I
forecasts are per-nuisance best-fit ranges; they are not complete allowed
envelopes or guaranteed discriminators. Parameter fits are not causal diagnoses.

The new numerical tests cover continuum gradients/covariance, Gaussian smoothing,
equivalent-width systematics, gap-safe integration/event screening, signal/nuisance
projection, physical units/geometry and normalization-invariant CO indices.
Native data checks and sensitivity experiments remain separate from these tests.
See each result subdirectory's report and protocol for the actual findings and
limitations, including inconclusive outcomes and rejected external assumptions.

The physical-model command also requires the preserved official SVO Cousins R
response at `models/filter_CousinsR.xml` in its result directory. Its exact
request and checksum accompany the data snapshot; it is not fetched silently
by the calculation. All commands above reuse existing native snapshots.
