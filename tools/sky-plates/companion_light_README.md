# Conditional companion-light shape comparison

Use `companion_light_sensitivity.py` from repository root with `python`. Add `--wide` for the broader primary-grid sensitivity. It consumes saved CK04, Gaia XP and SPHEREx inputs; no downloads. Dependencies are the existing sky-plates environment.

Results and predeclared protocol live in `results/sky-plates-j0754-companion-light-20261109`. Giant-only model adequacy is assessed before considering conditional composite-light envelopes. Every dataset is normalized/fitted separately, not stitched across epochs. Fraction is monochromatic 1 um light, not mass or bolometric luminosity. Log RMS tolerances are descriptive sensitivities, not likelihood/confidence bounds; native/local SPHEREx share photons. FAST is not included in quantitative continuum fitting because calibration is uncertified. Same-color components and spot-versus-companion interpretations remain degenerate; large-area dwarf-like templates need not be physical binaries. Saved grid excludes temperatures outside 3500-6500 K.

Run four numerical tests with `python -B -m unittest discover -s tools/sky-plates -p test_companion_light_sensitivity.py`. See the result report for numbers and limitations.


