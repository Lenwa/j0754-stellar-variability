# Existing TESS multi-day eclipse sensitivity

Run `python -B tools/sky-plates/multiday_eclipse_check.py` from repository root. This reads existing quality-screened native FITS, uses hourly averages and joint sector baselines, and tests 1–6-day periodic trapezoids without the earlier one-day highpass. No downloads or original-data writes.

Outputs and predeclared protocol are in `results/sky-plates-j0754-multiday-eclipses-20261109`. Native/local corrected products share photons. Phase/depth injection recovery is conditional on extracted products and cannot establish upstream attenuation or a universal eclipse/companion limit. Original-data structured-variation maximum is a descriptive recovery threshold, not a false-alarm probability. The report preserves phase gaps, nuisance sensitivity, period-error sensitivity and selected-trough non-recurrence.

Tests: `test_multiday_eclipse_check.py`; existing sky-plates dependencies/runtime. Source belongs to this existing registered repository; result data stay outside Git.
