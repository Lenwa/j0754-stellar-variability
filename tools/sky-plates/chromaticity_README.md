# Existing-data chromaticity checks

Run `python -B tools/sky-plates/check_existing_chromaticity.py` from the repository root. No network or original-data writes.

The script reuses the existing ASAS-SN parser, fits camera offsets/trends and independent V/g sinusoidal coefficients during their overlapping 2017/18 season, and records shared-time 7/14-day block-sign residual sensitivities. These ranges are descriptive, not calibrated confidence intervals. It also tests previously saved matched optical/W1 visits against endpoint-anchored affine-flux and linear-magnitude trajectories, excluding training years from residual summaries. No physical-model probabilities or robust hysteresis are inferred.

Inputs and protocol: outer results/sky-plates-j0754-chromaticity-20261109. Three tests verify exact coefficient recovery, gray-law equivalence and endpoint anchoring. Results remain outside source Git.
