# J0754 spectral-change comparison

This extends the existing historical sky-plate investigation. It compares
temperature/spot changes, uniform or partial dust obscuration, and two-component
occultation against the measured modern optical–infrared fading. It preserves
predictions for additional colors without selecting a cause in advance.

The maintained source is this directory in `this repository`.
Dependencies are the existing `requirements.txt`; use the isolated outer
`python`. Native observations, CK04 spectra,
SVO responses, protocols, model records and figures stay in outer
`results/sky-plates-j0754-models-20261108`, separate from Git source.

## Reproduction

From the repository root:

```powershell
$skyPython = 'python'
$skySource = 'tools/sky-plates'
$skyModels = 'results/sky-plates-j0754-models-20261108'
& $skyPython -B "$skySource/fetch_stellar_grids.py" --out "$skyModels/models"
& $skyPython -B "$skySource/compare_stellar_models.py" --inputs "$skyModels/models" --dust 'results/sky-plates-j0754-dust-20261108' --out $skyModels
& $skyPython -B "$skySource/compare_j0754_events.py" --historical 'results/sky-plates-j0754-20261107/dasch_differential.csv' --asassn 'results/sky-plates-j0754-independent-20261108/asassn_v1_lightcurve.csv' --legacy-asas "$skyModels/history/asas_annual_by_dataset.csv" --legacy-controls "$skyModels/history/controls/annual_exact_frame_differential.csv" --out $skyModels
& $skyPython -B "$skySource/orbit_photometry.py"
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

The legacy archive extension has its own [guide](historical_extension_README.md).
`orbit_photometry.py` uses the preserved earlier ASAS-SN/ZTF/TESS snapshots and
fixed output directory recorded in its source; it does not download data.

## Numerical contract and limitations

`stellar_models.py` integrates FLAM spectra through typed SVO energy/photon
responses, checks filter identity, and rejects incomplete wavelength coverage.
The generic SDSS g response approximates ASAS-SN g. A ZTF g sensitivity is
reported separately. Native CK04 grids are 250 K and relatively coarse in the
infrared; interpolation in 2/10 K steps creates no corresponding precision.
The baseline temperatures, gravities and metallicities are broad nuisance
grids, without a calibrated static-SED prior. Static foreground reddening is
not separately fitted into each baseline spectrum; its band-weighting effects
belong among the uncalibrated model systematics.

March 2018→2022 g=0.312 and W1=0.1175 mag are the primary training changes.
W2=0.08875 is excluded from optimization; its value was already known, so this
is a retrospective check. October values and response/atmosphere variations are
sensitivities. Two fitted quantities reproducing two training bands is not
evidence of a unique mechanism or of higher physical probability.

The JSON retains per-nuisance best fits separately from sampled compatibility
envelopes. Envelopes include all sampled opacity/temperature/area cases with
analytic covering/gray-offset intervals that keep both training bands within
±0.02 mag. This addresses the overly narrow ranges obtained by retaining only
one optimum per nuisance case. Neither range is a confidence interval, and no
W2 filtering is used to make the additional-color forecasts look stronger.

The conditional V−I and H−Ks predictions help design measurements, but their
smallest family gaps are comparable to the assumed per-band/model systematic
scale. Accurate matched-state measurements and multiple colors are needed;
they are not guaranteed discriminants. The broader g−Ks envelopes overlap.
Forecasts describe the historical bright/faint state change, not today's color.

Uniform foreground screens and purely gray single-component fading cannot fit
both training bands within the stated tolerance. Temperature plus scale, cool
spots, partial thick covering and two-component dilution each have viable
sampled cases. Fixed-area cooling has compatible cases at the coolest tested
baseline boundary; it is not excluded as a whole by the earlier blackbody toy.
Spot-temperature floor cases and large companion-area/covering values remain
explicit. A scale factor is not a radius measurement, and a two-component
spectrum is not a validated coeval binary or dynamical occultation model.

Eight new numerical tests cover detector integration, grid coverage, temperature
and geometric-scale recovery, physical partial covering, two-component dilution,
starting spot coverage, filter identity and non-optimal compatible envelopes.
The fixed-period test has separate cadence, nuisance-projection, coherence and
injection tests. A nonpersistent exact-period photometric signal does not
exclude companions, changing rotation signatures or other nearby periods.

## Seasonal activity-rhythm check

Run `check_activity_rhythm.py` from the repository root with the sky-plates
Python environment. It reuses the saved daily ASAS-SN g photometry and the
existing fixed-period fitting helper. Outputs are under
`results/sky-plates-j0754-activity-rhythm-20261108`.

This exploratory diagnostic imposes the Gaia 35.9204566-day period separately
in July–June seasons, jointly fitting camera offsets and camera-specific linear
trends. Seasonal amplitudes and phases vary; the 2024/25 fitted semi-amplitude
is only 0.0051 mag with a 0.0029 mag white-noise formal error. These fits do not
independently detect rotation. Correlated variability is not covered by the
formal errors, and phase is uninformative at small amplitude. The preserved
frozen-waveform validation (26% worse held-out squared error for the primary
ASAS-SN model) is the stronger test against a persistent photometric clock.

Same-camera annual medians retain the slow fading and recovery. Together with
the earlier spot/color feasibility calculation, this leaves evolving spot
coverage viable. It establishes neither a magnetic-cycle period nor companion
driving: one broad modern dip is insufficient, and brightness is not a direct
temperature measurement. The prior 912-day test also rejects that period as a
stable dominant explanation for the modern dip.

## Rotation search allowing seasonal changes

`search_seasonal_rotation.py` runs from repository root with the sky-plates Python.
It searches 10–100 days, allowing independent sine/cosine coefficients each
July–June season, with joint camera offsets and linear trends. Equal-season
mean fractional residual-SSE improvement peaks at 35.78 days in ASAS-SN (nine
seasons) and 36.30 days in ZTF (six qualifying seasons). Eight ASAS-SN seasonal
maxima lie at 34.1–37.3 days; the remaining weak season peaks at 40.1 days.
ZTF has five maxima at 35.4–37.3 days and a weak 2024/25 season at 41.9 days.
Thus a roughly 36-day recurring timescale survives changing seasonal phase.
The precise peak difference is not a measured differential-rotation rate.

Results, search protocol, per-season power curves and figure are in
`results/sky-plates-j0754-seasonal-rotation-20261108`. None of 256 block-sign
null draws at any of three block sizes (3/7/14 days) matches the observed
combined maximum in either survey; these are correlated-noise sensitivity
diagnostics, not calibrated false-alarm probabilities. An independent direct
joint least-squares calculation at three frequencies agrees with the ASAS-SN
2019/20 projected power to <1e-9. The search supports spotted rotation near
the candidate orbital period; it does not prove synchronization, confirm the
companion, identify the multi-year cycle, or model arbitrary within-season
spot evolution. No new data were downloaded.

## Rotational amplitude versus brightness

Run `compare_rotation_brightness.py` from repository root. It uses the same retained
seasons and daily data, comparing joint-camera fitted peak-to-peak amplitudes
with bq-only ASAS-SN medians and independent differential ZTF medians. Both
fixed 35.92-day fits and the previous search's selected 30–42-day seasonal
periods are retained. The latter amplitudes are subject to period-selection
bias. Results and figure: `results/sky-plates-j0754-amplitude-brightness-20261108`.

There is no consistent monotonic brightness–amplitude relationship: selected-
period descriptive Spearman coefficients are -0.033 (ASAS-SN, nine seasons)
and +0.086 (ZTF, six). These are not calibrated significance tests. Fixed-
period peak-to-peak amplitudes remain roughly 0.05–0.07 mag across several
bright/faint seasons, then fall to 0.010 (ASAS-SN) and 0.0048 mag (ZTF) in
2024/25. Nearby selected-period fits in that season reach about 0.022 mag;
do not confuse weakness of the exact 35.92-day component with no variability.
ASAS-SN's 2025/26 amplitude rebounds to 0.038 mag at the fixed period.

The simple larger-rotational-swing-when-fainter prediction is not supported.
Spots are still compatible: rotational amplitude measures longitudinal
asymmetry, not total coverage; interpreting the decoupling as a particular
spot distribution is an untested explanation. This adds evidence for changing
modulation, not a measurement of spot growth/decay or a magnetic cycle.
Direct fitted ASAS-SN amplitudes agree with the earlier independently projected
fixed-period fits to better than 1e-9 mag. The generated figure was inspected.

## Public-source attribution

### Shorter 2024/25 fits

`check_short_rotation.py`, run from repository root, fixes 35.9204566 days and uses
two calendar halves split on 2025 January 15 plus 90/120-day windows every
30 days from 2024 September 1. Per-camera offsets and linear trends remain
jointly fitted. Sparse or poorly conditioned fits are explicitly omitted.
Protocol and measurements: `results/sky-plates-j0754-short-rotation-20261108`.
Early 120-day amplitudes are 0.036 mag peak-to-peak (ASAS-SN) and 0.032 mag
(ZTF), above the corresponding whole-season 0.010 and 0.0048 mag. First-half
fits give 0.033 and 0.025 mag. Later ASAS-SN modulation weakens; ZTF suggests
phase changes but is sparse and its late amplitudes are less convincing under
the block-sign diagnostic. The season-long fit therefore understates some
within-season variation. Neither a persistent large hidden signal nor rapid
spot evolution is established: amplitude evolution, phase evolution and a
nearby period remain confounded. Overlapping windows are not independent;
512-draw seven-day sign tests are pointwise diagnostics, not trial-corrected
FAPs. A 0.05-mag peak-to-peak synthetic signal plus camera nuisance terms on
the actual early ASAS-SN cadence recovers its coefficients to <1e-9 mag.

Publication-context check (2026 November 8): AAS describes RNAAS as editor-
moderated and non-peer-reviewed, with 1,500 words and one figure or table.
https://journals.aas.org/research-notes/ . Exact-name web searches for the
target, Gaia ID and TIC ID with flare/TESS terms found no prior report of the
2021 December event; this is not a comprehensive ADS or flare-catalog event
cross-match and does not establish novelty. The archived TESS pixels are
public; our analysis has not been submitted or published. An outside scientific
review can assess a narrow flare/rotation report without resolving the long
dimming mechanism.

Atmospheres: [STScI CK04 atlas](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas).
Responses: [SVO Filter Profile Service](https://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=voservice).
Dust: [Draine grain models](https://www.astro.princeton.edu/~draine/dust/dustmix.html).
Native request receipts and hashes preserve the exact downloaded versions.
Gaia XP, LAMOST and SPHEREx audits are results, with measurement availability
distinguished from usable spectra or extracted target photometry.
