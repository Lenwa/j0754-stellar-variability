# J0754 flare comparison and evidence consolidation implementation plan

> For agentic workers: use superpowers:executing-plans for the authorized inline analysis, in the user's order: compare the four events, then assemble the evidence table.

Goal: measure the four saved TESS events consistently and reconcile all J0754 evidence into a source-linked table.

Architecture: preserve native inputs and earlier exploratory results. A new offline comparison script consumes six calibrated pixel cubes, the retained TESS response and the existing candidate windows. A separate curated evidence table links measurements to their native analyses and distinguishes independent observations from alternative reductions.

Tech stack: existing Python/NumPy/SciPy/Astropy/Matplotlib environment; no new dependencies.

Spec: user-approved five-step game plan, narrowed on 11 November to steps 2 then 1. The comparison covers duration, equivalent duration, conditional energy, localization and coverage/sensitivity. The evidence table includes successes, failed tests, limitations and remaining questions. Other-star screening remains the next separate stage.

Constraints: no changes to native telescope files; no forced FRED classification; no double-counted reductions; signed integrals; no interpolation across unobserved gaps without explicit model labeling; separate conditional geometry/temperature uncertainty from measurement noise; no email or submission. Existing public-publication authorization applies after privacy review.

Review focus: baseline curvature can mimic a broad event; unresolved sources are not distinguished by a source-centered centroid; sparse/flagged cadences cannot be counted as observed; empirical-template astrometry is not an official PRF posterior; old master report predates the additional candidates and eROSITA identification.

## Tasks

- [x] Add numerical safeguards before implementation: finite-exposure asymmetric profile integral; signed gap-aware area; duration and flux scaling; injected image-offset recovery; no live-time inflation or duplicate cadence coverage.
- [x] Implement `tools/sky-plates/compare_j0754_tess_events.py`. Fit finite-exposure asymmetric exponential rise/decay plus local baseline; compare apertures 1/1.5/2 pixels and linear/quadratic baselines. Report observed-sampled and full-model ED separately, FWHM and rise/decay, and the common blackbody conversion assumptions.
- [x] Calibrate empirical difference-template localization with nearby blank-time image injections, retaining failed fits and selection scope. Test coarse time-resolved displacement where signal permits; no blanket asteroid exclusion.
- [x] Reuse measured search injections and calculate a new direct-aperture sensitivity check at native retained timestamps. Report accepted cadence coverage and model-specific effective search time, not an unrestricted population flare rate.
- [x] Inspect figures and results. Reconcile any failure without changing the original results. Save input/source hashes and the analysis report under `results/sky-plates-j0754-flare-comparison-20261111`.
- [x] Only after the comparison, create a complete source-linked evidence table and short synthesis under `docs`, including all archive nulls and rejected explanations. Update repository entry points.
- [x] Run relevant numerical tests and the component suite, privacy scan and scoped diff checks. Publication uses only these scoped files under existing authorization; the final commit/remote state is verified after this document is committed and recorded in the delivery receipt.

Execution ruling: this is a bounded continuation of the existing research scripts in the clean, dedicated astronomy checkout. Preserve the registered repository and existing result paths; do not move the project or create another data copy. No parallel implementation is scheduled.

Review checkpoint:10new numerical safeguards and166existing tests pass(176total). Matched-control sampling was corrected after review; source/response/input hashes verify. Hourly astrometry remains unresolved for two events, while gradient-nuisance amplitude changes are below0.3% for event means. Baseline envelopes are prominent.36evidence entries link44local source files and37result families. Browser policy blocked file-URL preview; static HTML/image/source structure is verified without a browser workaround.

Completion checks:176/176tests passed; source/privacy scans clean;42HTMLtable rows,71links and2local images checked with no missing targets. No important review finding remains unaddressed. New result payload is approximately1.05MB. Scientific limitations remain in the evidence table rather than being treated as software failures.
