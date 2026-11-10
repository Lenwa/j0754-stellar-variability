# J0754 companion-light sensitivity, 9 November 2026

## Finding
A single giant atmosphere already fits the saved broad spectral shapes to approximately 1.74% log RMS in Gaia XP and 3.43-3.61% in each SPHEREx reduction/season. These are descriptive shape residuals, not uncertainty estimates. Adding a second atmosphere lowers residuals, but does not establish a companion-light detection. The allowed mixtures are strongly degenerate with temperature, reddening and calibration. These data do not yield a robust general upper limit on companion light.

## What was actually fitted
Existing local data only, without downloads. Gaia sampled XP was averaged over 20 nonoverlapping 30 nm bins across 400-1000 nm (actual sample means 400-998 nm), retaining its multi-epoch mean status. Model averages use exactly the same sample positions. This is an approximate broad-shape comparison, not a native XP response/covariance fit. SPHEREx uses the already declared primary-quality screen, then native flags=2^21 and fit_ql<5. Each 2025/2026 season is separate; native and local reductions are compared separately, never counted as independent photons. There are 35 selected points in 2025 and 33 in 2026, spanning 0.75-2.4 um. The primary fit excludes wavelengths beyond 2.4 um; the retrospective longer-wavelength diagnostic below evaluates 2.4-5 um without refitting.

Each dataset has its own free normalization, primary atmospheric parameters and extinction nuisance value. There is no combined multi-epoch likelihood or assumption of equal source brightness. Per-epoch fractional envelopes must not be intersected as if a constant companion fraction were expected: a constant-luminosity companion would change its fraction when the primary changes. Extinction variation here is a fitting nuisance, not evidence that foreground dust changes between epochs.

Models are saved CK04 atmospheres: initial primary 4000-5500 K in 250 K steps; log g  1.5,2,2.5,3; metallicity-1,-0.5,0. Companion 3500,4000,4500,5000,5500,6000,6500 K at log g  4.5 with the primary's metallicity. These are main-sequence-like spectral templates, not a self-consistent mass/radius/isochrone binary. Each spectrum is normalized at 1 um before mixing. The fitted companion fraction is fraction of monochromatic flux density at 1 um, not bolometric fraction, mass fraction or surface area. It is the same for F_lambda and F_nu at this wavelength. Grid fractions 0-.9 in.05 steps. Both components share DraineRv3.1 attenuation, Av 0-1.5 in.1 steps. SPHEREx responses are approximate rectangles, not native calibrated response curves. Spectra are mixed linearly before computing log residuals; normalization is the geometric mean observed/predicted flux ratio. The score is unweighted RMS log residual, not chi-square.

FAST was considered but excluded from quantitative broad-SED fitting: its absolute/relative continuum calibration is uncertified, native line-spread functions not known precisely, and the coarse atmosphere grids cannot turn its broad line comparisons into a precise composite-spectrum constraint. Existing Halpha changes retain their prior caveats; they are not converted into companion fractions. This is a limitation of the proposed three-source comparison, not evidence against a companion.

## Numerical results
Dataset                    Giant-only RMS       Mixture RMS (initial)
Gaia XP mean               .0173666              .0109457
SPHEREx native 2025         .0360692              .0317778
SPHEREx native 2026         .0356485              .0260455
SPHEREx local 2025          .0342863              .0287315
SPHEREx local 2026          .0346898              .0273604

Every giant-only model passes the exploratory 5% shape tolerance. At 3%, all SPHEREx giant-only fits fail; interpreting the composite models passing that arbitrary threshold as companion detections would mistake model inadequacy or calibration residuals for evidence. The known calibration/systematic issues do not permit a likelihood-ratio claim. The improvement uses extra freedom and was assessed on the same data being fitted, not held-out predictions.

At 5% tolerance, the largest sampled compatible 3500 K secondary fractions at 1 um are.50 for XP,.45 for native 2025,.40 for native 2026,.40 for local 2025 and.40 for local 2026. These are sampled conditional compatibility endpoints, not upper confidence limits. At 10% they become.70,.90,.85,.90,.85; .90 is the grid edge, not an exclusion above it. Fractions for 4500-5500 K secondaries reach.90 for every dataset at 5%. For 6500 K they reach.85 in XP and.90 in SPHEREx. Consequently there is no useful general small-light-fraction conclusion. Do not summarize these as 'the companion contributes 40%' or a measurement of its temperature.

## Boundary and degeneracy checks
The first SPHEREx composite optima put the primary at the 4000 K boundary. A documented post-initial-run amendment expanded primary temperatures to all locally saved 3500-6500 K. Giant-only results remained unchanged; composite best RMS became.0315656,.0245925,.0286082,.0254375 for the four SPHEREx datasets. Some composite fits still touch cool/hot boundaries. The key 3500 K compatibility endpoints at 5% were unchanged; other broad intervals expand. Exact expanded profiles are in results_wide.json.

An identical spectral component is mathematically indistinguishable from a rescaling, at any mixture fraction. The unit test confirms this exactly. Similar-color components also hide readily. Allowing an arbitrary emitting area means some compatible large-fraction log g  4.5 templates could be physically implausible stars; the shape-only experiment deliberately does not promote them to coeval binaries. A companion mass-radius relation, primary radius/distance constraints, realistic atmospheres and consistent photometric states would be needed to restrict those cases.

The composite mathematical form can also describe two temperatures on one spotted star. Different gravities/line profiles can distinguish those physical situations with sufficient calibrated spectral information, but broad-shape residual improvement alone cannot. This study neither establishes a luminous companion nor weakens the existing activity/spot evidence on its own. A faint, dynamically relevant companion remains fully compatible. Companions hotter than 6500 K, cooler than 3500 K or compact-remnant atmospheres were not tested and are not excluded.

## Saved artifacts and reproduction
Maintained code: tools/sky-plates/companion_light_sensitivity.py and test_companion_light_sensitivity.py in build/github-desktop-review.
Run from outer Lumen: runtime/sky-plates/venv/Scripts/python.exe -B tools/sky-plates/companion_light_sensitivity.py ; add --wide for primary-grid sensitivity.
Results: results.json, results_wide.json, shape_comparison.png, shape_comparison_wide.png, input_hashes.json, protocol.txt. Input hashes cover all saved atmosphere files, Gaia CSV, both SPHEREx input records and dust table. No original data were modified. Figures visually inspected. Four numerical tests passed: scale invariance, known mixture recovery, identical-template degeneracy and positive-flux guard.

No new observations, contact, cloud/API charges or publication. No push. Source commit status is reported separately; native data/results remain outside source Git.

## Retrospective longer-wavelength check
After the primary comparison, the selected models were evaluated against 2.4-5 um SPHEREx points with their original nearIR-derived normalization, without refitting or excluding any compatible models. These data were already known, so this is retrospective, not a blind prediction. For wider-primary results, median observed/model ratios in 3.8-4.4 um are:
Native 2025: giant 1.0312, mixture.9614; local 2025: giant 1.0100, mixture.9287.
Native 2026: giant 1.0297, mixture 1.0037; local 2026: giant 1.0306, mixture.9953.
In 4.4-5 um the 2026 values are giant 1.1011 versus mixture 1.0582(native), and 1.0930 versus 1.0374(local). Thus mixture improvement in nearIR does not consistently improve longer wavelengths across epochs/reductions. Known atmosphere/CO/response/calibration limitations prevent assigning these residuals to companion light or dust. Full band counts, ratios and fixed-scale RMS are in longwave.json/longwave_wide.json. No saved compatibility envelope was tightened by these checks.

Plot caption: Gaia panel displays normalized F_lambda; SPHEREx panels normalized F_nu. Each panel has its own scale. Similar plotted shapes are not a common absolute flux scale.

## Explicit size illustration, not a fit
Using the saved solar-metallicity CK04 surface flux at 1 um, an ASSUMED 3500 K/log g  4.5 dwarf of radius.3 Rsun next to an ASSUMED 9 Rsun giant(log g  2.5) has companion/primary flux ratio.000525 for 4500 K primary,.000361 for 5000 K and.000309 for 5250 K. Total-light fractions are essentially.0525%,.0361%,.0309%. This illustrates why a small companion can remain invisible in broad spectra even if dynamically relevant. Neither radius is measured here; these values are not inferred from Gaia mass function and are not a luminosity bound or prediction for the actual companion. The faint lowmass hypothetical companion could also be cooler than this grid's 3500 K floor.


