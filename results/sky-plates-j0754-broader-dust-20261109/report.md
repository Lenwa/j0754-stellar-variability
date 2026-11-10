# J0754 broader dust tests using existing observations

9 November 2026. All calculations use existing local observations, CK04 spectra, Draine extinction tables and saved filter responses. No new observations or downloads were obtained.

## Main result

Broader dust explanations can be tested, and several survive. The existing three-band fading does not securely distinguish cool spots from partial-covering or sufficiently flat extinction. One new temporal test disfavors a specific simple dust history: a fixed-area cloud that starts transparent and changes only its column predicts the intermediate infrared visits substantially worse than changing coverage or fixed-temperature spots. A heavily obscured baseline improves that temporal fit, but the specific tested rescue is much too red in absolute g-W1 under its single-giant assumptions. Other dust geometries remain viable.

## What was already tested

The earlier integrated model comparison already included partial covering for Draine Rv 3.1, 4.0 and 5.5, not only uniform dust screens. All 108 atmosphere/law nuisance cases fitted g and W1. Their best-fit W2 predictions were 0.0703-0.0750 mag, compared with 0.08875 observed. Residuals of -0.0185 to -0.0137 mag survive the adopted descriptive systematic allowance. The prior full two-band 0.02-mag tolerance envelope extended to W2 0.0566-0.0899 mag. These were already viable alternatives.

Ordinary uniform screens failed the g/W1 comparison (0 of 108 compatible cases; best maximum training residual about 0.1005 mag). The prior illustrative spot prediction was W2 approximately 0.0995 mag, residual +0.0108. That modest difference from partial dust is not a reliable model ranking. The broader spot grid had W2 residuals from approximately -0.0305 to +0.0375 mag among training-compatible nuisance fits.

## Method and new model families

The primary anchors are March 2018-to-2022 fades g=0.312 and W1=0.1175 mag. W2=0.08875 mag is excluded from parameter optimization and inspected afterward. It was already known, so this is a retrospective withheld-band check. October anchors (0.3225, 0.115, 0.090 mag) are a sampling sensitivity.

Attenuation is integrated through g, W1, W2, V, r, I, J and Ks responses with CK04 spectra. Standard-law calculations use 4500, 5000 and 5500 K, log g=2.5 and metallicity -0.5. The smooth power-law extension uses the nominal 5000 K spectrum. No thermal or scattered dust emission is added.

For partial covering, F_j(t)/F_intrinsic,j = 1-c(t) D_j(tau(t)), where j denotes band, t denotes epoch, and D is the band-integrated blocked fraction. F_faint/F_bright is computed explicitly, including nonzero initial coverage/column. We tested:

- Fixed optical depth, changing coverage, with initial coverage 0.0, 0.1, 0.3, 0.5 and 0.7, and tau at 0.48 microns sampled from 0.0001 to 300.
- Fixed coverage, changing column, with initial tau 0.0, 0.1, 0.3, 1, 3 and 10. For each final tau, coverage is solved analytically from g: c=(1-r_g)/(D_faint,g-r_g D_bright,g). This avoids missing narrow solutions near the observed 25% optical loss with a coarse coverage grid.
- A gray attenuation term plus a standard selective extinction screen. The two training bands determine the two parameters.
- A phenomenological uniform power-law screen tau(lambda) proportional to lambda^(-alpha), with alpha between 0 and 2.5, and a partial-cover extension with initial coverage 0.0 or 0.3. A successful slope is not a measured grain-size distribution.
- Increasing cool-spot coverage plus a uniform Rv 3.1 foreground screen. Hot and cool atmospheres are each attenuated before integration. Spot temperatures span 3500 K to below the primary temperature in 250 K steps, initial spot fractions 0.0, 0.1 and 0.3. Final coverage is solved from g for each dust column; W1 then tests the remaining parameter.

For grid families, g is anchored exactly and the sampled W1 tolerance is +/-0.02 mag. These are conditional slices, not full two-band tolerance volumes or confidence regions. Gray-plus-selective and uniform power-law results below are exact two-anchor fits, not equivalently explored tolerance envelopes. Parameter grids and compatible samples are retained in results.json.

## Withheld W2 results

Observed W2 fade is 0.08875 mag.

| Family | Prediction from g/W1 | Interpretation |
|---|---|---|
| Uniform flat power law | 0.09998 mag | Residual +0.01123 mag, comparable to the illustrative spot model |
| Gray plus standard screen | 0.11287-0.11324 mag | Residual +0.02412 to +0.02449 mag; worse than the illustrative spot fit, but sensitive to the adopted systematic allowance |
| Nonzero-baseline fixed-tau varying coverage | 0.05049-0.08798 mag | Sampled exact-g/W1-tolerance envelope; compatible edge approaches observed W2 |
| Fixed coverage with changing column | 0.05374-0.08798 mag | Similar conditional envelope; temporal behavior supplies an additional test below |
| Partial power-law cover | 0.04651-0.12002 mag | Broad compatibility, not a unique extinction law |
| Spots plus standard screen | 0.06863-0.12983 mag | Mixed models remain underdetermined |

The uniform power law has alpha=0.50575 and tau(0.48 micron)=0.28692. It predicts V=0.29022, r=0.27419, I=0.24296, J=0.19306 and Ks=0.14581 mag for the same state change. This is much flatter than the tested ordinary interstellar laws across the optical-to-IR interval. Whether a physical grain population and geometry can realize it must be tested separately.

An illustrative mixed model has a 5000 K photosphere, 3500 K spots, spot coverage increasing from 0.3000 to 0.46464, and added tau=0.03938. It predicts g=0.31200, W1=0.11768 and W2=0.09261 mag. The small W2 residual (+0.00386 mag) demonstrates compatibility, not evidence that both mechanisms are present. Different spot temperatures, initial coverage and dust columns trade off.

The October sensitivity gives uniform power-law W2=0.09697 versus 0.090 observed, and gray-plus-standard W2=0.11006-0.11045. The qualitative conclusions are unchanged.

## New temporal test of a changing-column cloud

We reused the saved per-camera 2018 and 2022 g/W1 anchors and predicted all other matched W1 visits. Each model fits only the two anchor contrasts. With coverage and initial column held fixed, each evaluation visit's column is inferred from g alone; no evaluation W1 measurement chooses the fit. Camera streams share IR measurements and are not independent samples. Matches span +/-30 days, with the existing +/-15-day sensitivity retained.

For the best-sampled bq camera (10 evaluation visits), a cloud initially transparent (tau0=0) predicts W1 with RMS residual 0.0406-0.0428 mag across the three standard laws. Earlier fixed-temperature spot/fixed-opacity variable-covering trajectories have RMS 0.0206 mag, and a free-slope magnitude-linear relation has RMS 0.0191 mag. The initially transparent changing-column model is therefore a poorer description of these visits. Other cameras show the same direction: transparent-cloud RMS 0.0394-0.0427 mag in bi and 0.0349-0.0367 mag in bm, versus affine-flux RMS 0.0198 and 0.0208 mag.

Allowing tau0=3 improves bq residuals to 0.0165-0.0186 mag, but only 15.5-15.9% of intrinsic g flux reaches the observer at the bright baseline in those solutions. The intrinsic g luminosity would be about 6.3-6.5 times the observed baseline value. Initial tau is scanned as a nuisance assumption; choosing the lowest evaluation RMS afterward would be retrospective model selection, not independent validation. This is an additional substantial SED/luminosity/geometry requirement, not an inferred bolometric correction or a validated cloud. The additional baseline-color check below tests this requirement directly; the full noncontemporaneous SED has not been jointly fitted.

All 64 temporal model records retain every available evaluation visit; no visits fell outside their optical domains. All model evaluation rows and any out-of-domain optical visits are explicitly retained. There is no extrapolation to negative dust optical depth or silent clipping. The finite grid and matched-visit systematics prevent turning RMS differences into formal odds or a universal dust rejection.

## Conditional baseline color check

A targeted retrospective check tests the heavily obscured rescue against the already matched 2018 absolute color, g-W1=12.621-9.50675=3.11425 mag. We convert integrated CK04 F_lambda to equivalent band F_nu by dividing by the same passband integral of c/lambda^2. We use nominal AB g=3631 Jy and the saved SVO W1 Vega zero point=309.54 Jy as a nominal Vega flux-density proxy, not exact synthetic photometry against a calibrated Vega spectrum. Response/color corrections are not independently calibrated. The SVO g file defaults to Vega (4023.573 Jy); that value was explicitly not mistaken for an AB zero point. The ASAS-SN public documentation states calibration to APASS; the SDSS g response remains a proxy. This is a large-discrepancy diagnostic rather than precision synthetic-photometry validation.

The same 5000 K, tau0=3 temporal-rescue geometries predict baseline g-W1=4.4768-4.4857 mag, or 1.3626-1.3714 mag redder than observed. Keeping those geometries and changing the atmosphere to 4500 or 5500 K gives discrepancies of approximately 2.05 and 0.81 mag. These temperature variations are sensitivity checks, not refits. Adopting the saved g Vega zero point would shift predictions another +0.1114 mag, not remove the mismatch. Additional ordinary foreground reddening also moves in the wrong direction.

Thus the specific heavy-baseline single-giant rescue fits the relative trajectory but fails this absolute-color consistency check under its stated spectral assumptions. It should not be called physically established merely because its temporal RMS is low. A different intrinsic spectrum, additive/scattered light, or more complex geometry would be another model requiring its own evidence. This check does not reject all partial-covering models or the flatter uniform-screen alternative. Source: https://www.astronomy.ohio-state.edu/asassn/public/examples.shtml and saved SVO response metadata; numerical results in baseline_color.json.
## What the tests identify

Fixed-opacity variable covering and a fixed two-temperature surface both produce affine flux trajectories. Their identical functional form prevents the current optical/W1 track alone from distinguishing them. Changing column is different and is constrained by the new test. A substantially obscured baseline restores temporal flexibility, but the tested heavy-baseline rescue fails the absolute-color consistency check.

Two training bands determine at most two independent continuous combinations. Standard partial dust fits use two parameters for fixed baseline choices; power-law partial covering and mixed dust/spot interpretations introduce additional nuisance freedom. Three measured fades cannot determine geometry, extinction law, baseline obscuration, temperatures and spot fractions simultaneously. A low residual after adding freedom is not model evidence.

The ordinary uniform-screen failure remains strong. The broader statement that dust cannot explain the measured optical/IR changes is not supported. Spots remain physically motivated by the activity diagnostics, but partial covering, flat extinction and mixed models pass meaningful existing-data checks. New simultaneous optical colors and temperature-sensitive spectroscopy across brightness states would distinguish these surviving models more effectively than another two-endpoint fit.

## Reproducibility

Maintained source: broader_dust.py and test_broader_dust.py in the existing sky-plates component. Results retain every fitted nuisance case, conditional compatible sample, temporal prediction and source hash. No original data or report was edited. Software checks verify numerical implementation, not a physical origin.




