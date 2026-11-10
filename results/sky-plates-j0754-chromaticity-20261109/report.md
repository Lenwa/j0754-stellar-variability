# What existing colors can determine about J0754

9 November 2026. Existing-data tests only; source photometry unchanged.

## Result

The existing overlap shows the approximately 36-day modulation in both V and g with similar amplitudes, but does not detect a reliable difference between the bands. The full optical/infrared trajectory is roughly consistent with the known color-dependent fading, with deviations from a fixed two-component trajectory at the few-hundredths-of-a-magnitude level. Those deviations do not distinguish spots from other mixtures or establish hysteresis.

## Inventory and genuine overlap

The saved raw ZTF target file has 575 g measurements (470 accepted by the existing quality function), zero r and zero i. It cannot provide a g-r curve.

ASAS-SN has 960 accepted V exposures from 2012-02-16 to 2018-11-29 and 4220 accepted g exposures from 2017-10-07 to 2026-11-08. Nearest opposite-band matches within half a day exist for 85 V exposures on only 30 distinct V dates; just six exposure pairs are within 0.05 day. These counts do not mean 85 independent simultaneous colors. There is no V coverage of the 2022 minimum.

## New common-season V/g test

The protocol selected 2017-10-07 through 2018-05-01 before fitting. Accepted observations were reduced to camera-day median magnitudes, with no brightness clipping. Cameras needed at least ten daily groups. Retained data: 24 ba and 13 be V camera-days; 59 bq, 79 bi and 27 bm g camera-days. Camera-days across cameras can share a date and are not independent observing nights.

Joint equal-weight least squares uses a separate offset and linear trend for each camera, with separate V and g sine/cosine coefficients at the fixed 35.9204566-day period. The design is full rank (14 columns). A known injected coefficient vector is recovered to better than 1e-15 mag.

- V semi-amplitude: 0.03379 mag.
- g semi-amplitude: 0.03150 mag.
- g minus V cosine/sine coefficients: -0.00075, +0.00426 mag.

Seven-day shared-time block-sign residual perturbations give 2.5th–97.5th percentile amplitude ranges V 0.0252–0.0443, g 0.0288–0.0346 mag. The corresponding color-coefficient ranges are [-0.00956,+0.00885] and [-0.00319,+0.01180] mag; both encompass zero. Fourteen-day blocks give the same qualitative conclusion. These are descriptive correlated-residual sensitivities, not calibrated confidence intervals or a formal test of a zero color signal.

Consequently we cannot distinguish a modest temperature-dependent optical amplitude from gray modulation using this overlap. Similar V/g amplitudes do not prove ellipsoidal geometry: these nearby bands can respond similarly to cool spots, and the sparse V data, camera calibration and evolution within a season limit this test. It does add a second-band consistency check on the recurring modulation.

## New full optical/W1 trajectory test

Previous owner work already measured the 2018-to-2022 g/W1/W2 changes and same-frame infrared color. We did not claim those endpoint measurements as new. Here we asked whether the intervening and recovery visits follow a fixed spectral mixture.

For each camera separately, median matched 2018 and 2022 values anchor two descriptive trajectories: (1) affine relative flux, appropriate to changing weights of two fixed spectra, and (2) a linear magnitude relation. Neither family is a uniquely identified physical mechanism. Fixed-temperature spots and some partial-covering geometries share affine-flux behavior. The magnitude-linear law is free in slope; it is not a validated ordinary interstellar extinction law.

Only visits outside 2018 and 2022 contribute evaluation residuals. These data and endpoint years were already known, so this is a retrospective comparison, not blind validation. The existing plus/minus 30-day optical matches yield 7, 5 and 10 evaluation visits for bi, bm and bq. Affine-flux W1 residual RMS is 0.0198, 0.0208 and 0.0206 mag; magnitude-linear RMS is 0.0183, 0.0195 and 0.0191 mag. Their predictions differ by at most 0.0022–0.0028 mag at those visits, substantially less than the residuals. We cannot select one using these data. Camera comparisons share the same infrared measurements and are not independent replications.

For bq, the largest affine residual is -0.0392 mag in October 2020: W1 is brighter than the endpoint-anchored relation predicts. October 2020 and October 2023 have nearly equal matched g values (12.8110 and 12.8055) but W1 differs by 0.0305 mag. That is an interesting trajectory discrepancy, not a verified hysteresis measurement. The optical medians span weeks rather than matching each infrared exposure, the roughly 36-day rotational phase is unevenly sampled, and infrared/reference/camera systematics are not fully calibrated.

The tighter plus/minus 15-day selection leaves only bi with both anchor years and three evaluation visits; its affine RMS is 0.0131 mag. The better-sampled bq series loses the 2018 anchor, so it cannot supply the same fully reanchored tighter-window test. This is a real data limit, not a reason to manufacture a missing baseline.

## What follows

We can establish that the short modulation appears in two existing optical bands during one overlap season, and that the multiyear color trajectory is not a perfect fixed two-spectrum track at few-hundredths-mag precision. We cannot turn either finding into an orbital identification, a companion-light measurement, a spot-temperature measurement or an exclusive choice of fading mechanism. The target's non-gray long dip remains established by earlier work; this follow-up does not strengthen ordinary uniform dust into a competitive explanation.

## Saved state

Protocol, input hashes, coefficients, block sensitivities and all trajectory residuals are retained in this folder. Maintained source and three numerical tests are in the existing sky-plates component as check_existing_chromaticity.py and test_existing_chromaticity.py, with chromaticity_README.md. No downloads, contact, accounts, observing requests or source-data changes occurred.
