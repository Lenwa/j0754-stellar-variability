# FAST relative-velocity information audit

Run `fast_relative_rv.py` with the existing outer `python` from any directory. It locates the repository root root from the maintained source path, reads ten saved extracted spectra under `results/sky-plates-j0754-spherex-20261108/literature/fast`, and writes diagnostic results to `results/sky-plates-j0754-fast-rv-20261109`.

This is a relative blended-line-shift/consistency check, not a calibrated radial-velocity or orbit fitter. It keeps grating groups separate, records native hashes and BCV metadata, compares disjoint windows and same-night repeats, and tests continuum/template/smoothing sensitivity. BCV-added shifts assume an observed-frame native wavelength axis and remain conditional. Sky and telluric-region diagnostics do not certify stellar slit-centering corrections. Original spectra are never written.

Tests: `python -m unittest discover -s tools/sky-plates -p test_fast_relative_rv.py` using the existing astronomy environment. Dependency set is the component's existing numpy/scipy/astropy requirements. The report and protocol are research results outside source Git. The code's numerical injection success is not evidence of real RV accuracy.
