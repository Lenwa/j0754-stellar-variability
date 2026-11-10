# Existing TESS oscillation information

Run `python -B tools/sky-plates/check_tess_oscillations.py` from the repository root. The script locates that root from its maintained source path, reads the eight saved J0754 TESS FITS, and writes diagnostics to `results/sky-plates-j0754-oscillation-check-20261109`.

This is an exploratory relative-power comparison, not an asteroseismic pipeline. SAP and corrected sectors remain separate; the known flare is excluded and tested separately. Deterministic mode-envelope injections check final linear detrending transfer, not mission-pipeline or stochastic-oscillation detectability. No nu_max, seismic radius, suppressed-oscillation classification or upper limit should be inferred from absent visual features.

Tests: `python -B -m unittest discover -s tools/sky-plates -p test_tess_oscillation_information.py`. Existing numpy, scipy, astropy and matplotlib dependencies are reused. Native data/results remain outside source Git.
