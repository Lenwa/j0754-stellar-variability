# J0754 rotation cadence/alias diagnostic

Run from the repository root with python -B tools/sky-plates/check_rotation_aliases.py. The unchanged inputs are the prior daily ASAS-SN table and corrected ZTF measurements. Results and predeclared protocol are in results/sky-plates-j0754-alias-check-20261109, outside source Git.

The diagnostic compares the exact Gaia candidate period with its plus/minus annual-frequency aliases, half/double periods and the synodic lunar period. It preserves the previous seasonal selection, weighted camera offsets/slopes and independent seasonal sine/cosine coefficients. Noiseless phase-varied sinusoid injections on actual sampling measure cross-talk; they are not astrophysical false-alarm simulations. No new observations/downloads, confidence intervals, differential-rotation measurements or orbit confirmation follow. The previous broader search remains the source for the independently selected photometric peaks.

Three regression tests cover true-period recovery with nuisance terms, additive-offset invariance and annual-frequency arithmetic. Scope and outcomes are in the local report.
