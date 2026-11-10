# J0754: information in the existing TESS oscillation spectrum

9 November 2026. Existing-data diagnostic, not a seismic detection or nondetection limit.

## Finding

The saved TESS products have sufficient cadence to examine hours-long variability, but this bounded analysis finds no localized oscillation-like power envelope reproducible across sectors and reductions. It does not establish a seismic radius, mass, evolutionary state, or suppressed-oscillation classification. In particular, it does not justify applying the binary fractions for nonoscillating active giants from population studies to J0754.

There is a useful concrete warning in the actual data: leaving the already identified December 2021 flare-like event in sector 46 increases integrated 30–100 microHz power by **1.88 times in SAP and 3.04 times in corrected QLP flux**. A broad power excess in that sector can therefore be transient contamination rather than oscillations.

## Method and available information

All eight previously saved FITS were read, producing 16 SAP/corrected product comparisons. QUALITY=0 and finite positive flux were required. No brightness clipping was used. Each product/sector was normalized by its own median; a constant and linear trend were projected out. The previously localized flare interval BTJD 2571.5–2573.0 was excluded before computing the main spectrum, then retained as an explicit sector-46 sensitivity. No flux baselines were stitched.

Lomb–Scargle power was evaluated from 10 to 200 microHz in 0.05-microHz steps. For descriptive plotting, periodogram power was scaled by twice the effective mean cadence and converted from fractional flux squared per Hz into ppm squared per microHz (net factor 1e6), then smoothed with a Gaussian of sigma 5 microHz. This is an approximate PSD-like relative-shape diagnostic for irregularly sampled data, not a calibrated noise likelihood. A regular-cadence test verifies integrated sinusoid variance in these units. The figure and exact band powers are retained.

The observed native cadences are about 10 minutes in QLP sectors 44–46 and 3.33 minutes in sectors 71/72/1751. Typical long sectors supply about 2,700–9,000 retained points per product. The special sector 1751 has only 1,502 retained points and a much shorter baseline. The relevant giant oscillation timescales are below these Nyquist limits.

Illustrative solar scaling, nu_max = 3090 microHz times (M/R^2) divided by sqrt(T/5772), puts the assumed M=1–2.6, R=9.1–12.2 and T=4500–5206 range at roughly 20–110 microHz. These inputs are not a measured stellar solution; the scaling only defines a broad region to inspect. No measured nu_max or mode spacing was extracted.

## Observed power shapes

Most maxima of the smoothed spectrum within 25–110 microHz occur at the lower edge, consistent with a declining low-frequency background rather than a localized envelope. Corrected sector 46 has a descriptive maximum near 51.3 microHz and corrected QLP sector 71 near 67.05 microHz. These are not stable shared peaks across the other sectors or the alternate sector-71 SPOC reduction. Selecting one of them as nu_max would be unwarranted.

Native SAP low-frequency power differs strongly from corrected products. QLP and SPOC reductions of the same 2023 photons also differ. These comparisons measure processing sensitivity, not independent observations of different stars. Unknown granulation, activity, instrumental correlations and dilution are not separated by this calculation.

Typical robust adjacent-difference scatter is about 1,220–1,280 ppm in the 2021 QLP products, 2,290–2,450 ppm in 2023 QLP, and about 1,690–2,090 ppm in most SPOC products. This point-to-point statistic is not an oscillation-amplitude upper limit and does not represent red noise.

## Synthetic signal transfer

A deterministic mode-comb illustration uses frequencies 24–96 microHz spaced by 6 microHz, Gaussian mode amplitudes centered at 60 microHz with sigma 15 microHz, and seeded random phases. Its RMS is set to 100, 300 or 1,000 ppm at each actual retained cadence. Linear detrending preserves essentially all of its RMS, so the absence of a clear hump is not explained by this particular final detrending step removing such signals.

For corrected products, the synthetic 100-ppm signal alone contributes roughly 2–6% of the observed integrated 30–100 microHz power. At 300 ppm that ratio is roughly 19–55%. These are signal-power leverage comparisons, not detection probabilities: a coherent comb may be measurable below total broadband noise, while real solar-like modes have stochastic lifetimes and amplitudes. The test injects after upstream mission reductions and cannot measure their attenuation. It therefore cannot establish how much real oscillation power would have been recoverable, or label the star nonoscillating.

## What this determines

- The existing cadence contains relevant frequency information; lack of suitable cadence is not the obstacle.
- The known transient materially contaminates that frequency range.
- Processing-dependent background dominates this simple diagnostic; no consistent envelope is identified.
- No new seismic constraint follows from this bounded calculation. It neither confirms nor rejects the companion/activity explanation.

## Reproduction

Maintained source: tools/sky-plates/check_tess_oscillations.py and test_tess_oscillation_information.py. Run the source with the existing runtime/sky-plates/venv/Scripts/python.exe. Results, original-input hashes, protocol, power arrays and figure are in this folder. Original inputs unchanged; no downloads, contacts, observations or publication.

Four numerical tests cover trend projection, mode-envelope transfer, frequency recovery and integrated power scaling. Final suite and source commit status are reported in verification.json.
