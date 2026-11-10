# J0754: fixed 35.9204566-day photometric test

The tested stable waveform does not persist across independent epochs. Applying the training waveform worsens held-out ASAS-SN and ZTF fits; TESS cannot establish a coherent counterpart from these reductions and sector baselines.

The prediction/protocol was saved before fitting. No period scan or brightness clipping was performed. No external data were downloaded.

## ASAS-SN g: primary bounded test

5180/5249 valid individual measurements across V+g; 1396 daily camera groups in g.

|Joint nuisance|Full semi-amplitude (mag)|Train amplitude|Validation amplitude|Phase difference (cycles)|Held-out SSE improvement|
|---|---:|---:|---:|---:|---:|
|spline365|0.01136|0.02566|0.00068|-0.289|-26.02%|
|spline180|0.01073|0.02463|0.00087|-0.412|-25.97%|
|spline730|0.01064|0.02464|0.00061|-0.242|-19.42%|
|season|0.01064|0.02473|0.00105|-0.381|-26.78%|
|equal_weight_spline365|0.01407|0.02770|0.00419|-0.140|-8.77%|
|spline365_plus_harmonic|0.01134|0.02566|0.00071|-0.286|-26.18%|

Training ends 31 December 2021; validation begins 1 January 2022. Phase is the cosine-maximum convention at JD 2457000. Each validation fit re-estimates nuisance terms; held-out coefficients remain frozen to training. Positive held-out improvement alone is insufficient for a detection. Validation amplitude is consistent with zero, so its point-estimate phase is poorly defined.

Primary full-series amplitude 0.01136 ± 0.00153 mag (white-noise-only 1σ); validation 0.00068 ± 0.00186 mag. Train/validation coefficient-discrepancy χ²=79.10 for two coefficients; this is not calibrated for red noise.

## Actual-cadence injection/sensitivity

The primary nuisance model retains 0.990–0.994 of signal energy over all sine phases. Correct joint fitting recovers paired coefficient increments to numerical precision, while arbitrary detrending before fitting would attenuate the signal.

200 within-camera circular shifts of the real null residual yield median/95th-percentile apparent semi-amplitude [0.00268957677013902, 0.0051346085272784454]. These shifts preserve camera sample order but distort irregular time spacings; they are a sensitivity diagnostic, not a false-alarm calculation.

|Injected semi-amplitude|Recovered total amplitude range over 12 phases|Maximum total phase error (cycles)|
|---|---|---:|
|0.010 mag|0.00139–0.02136 mag|0.468|
|0.030 mag|0.01864–0.04136 mag|0.061|
|0.100 mag|0.08864–0.11136 mag|0.018|

## Independent ZTF g diagnostic

- spline365: full semi-amplitude 0.01840 mag; held-out SSE improvement -63.19%; phase difference -0.232 cycles.
- season: full semi-amplitude 0.01760 mag; held-out SSE improvement -64.59%; phase difference -0.244 cycles.

ZTF uses existing same-exposure differential estimates and equal daily weights; no per-point uncertainty was present in the retained JSON. V-band is a separate diagnostic in structured results, split at 2017 because most V coverage precedes 2022. ASAS-SN V1 has no explicit image-quality flags; the retained parser rejects limits, invalid errors/flux, sentinels and inconsistent camera filters. There is no brightness clipping.

## TESS: short-sector ambiguity

QUALITY==0, finite positive flux. All eight files were inspected. QLP SAP and corrected KSPSAP/DET flux, SPOC SAP and PDCSAP flux remain separate. Some are duplicate reductions of the same 2023 epochs. Per-sector medians and fitted offsets/slopes remove longterm information. The fixed 36-day signal has less than one cycle per usual sector; short-sector inferred amplitudes can be strongly degenerate with trends. Joint scores use effectively equal daily weights (0.02 mag floor); individual sector diagnostics use a 0.001 mag floor. Sector 1751 has only five valid daily groups from 18–22 January 2026 and is insufficient for its own periodic fit.

|Reduction|Sectors|Full semi-amplitude|Phase-energy retained range|Independent 2021 vs 2023+ comparison|
|---|---|---:|---|---|
|QLP_SAP_FLUX|44,45,46,71,72|0.00144 mag|0.225–0.348|held-out improvement -45.5%, phase shift +0.439|
|QLP_KSPSAP_FLUX|44,45,46|0.00013 mag|0.208–0.408|insufficient independent early epoch|
|QLP_DET_FLUX|71,72|0.00021 mag|0.221–0.292|insufficient independent early epoch|
|SPOC_SAP_FLUX|1751,71,72|0.02358 mag|0.200–0.237|insufficient independent early epoch|
|SPOC_PDCSAP_FLUX|1751,71,72|0.00309 mag|0.200–0.237|insufficient independent early epoch|

## Interpretation limits

A stable signal at the exact Gaia candidate period is not established when epoch phase and held-out benefit disagree or amplitude changes materially with nuisance/reduction choice. All quoted covariance uses independent daily residuals and excludes correlated/systematic noise, calibration drift, and Gaia period/ephemeris uncertainty. The fixed-period failure would not exclude a nearby drifting signal, rotation, an eclipsing system with non-sinusoidal or intermittent variability, or a companion without a photometric signature. Those hypotheses were not searched here.

The candidate Gaia orbit fails recommended gold-sample significance/efficiency cuts; even a photometric match would not establish a companion or the cause of the historical dimming. Maintained source is `tools/sky-plates/orbit_photometry.py`; results are local and uncommitted. No online publication.
