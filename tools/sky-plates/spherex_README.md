# SPHEREx extraction and archive extension

This component extracts J0754's public SPHEREx recovery-era spectrum and checks
its broad shape against stellar atmospheres. It extends the existing sky-plate
investigation; it does not identify the cause of fading. Maintained source is
`tools/sky-plates` in this repository. Native data, protocols, calibration maps,
literature and figures belong in the outer `results` directory, not Git.

Use the existing `runtime/sky-plates/venv` Python environment, with pinned
dependencies in `requirements.txt` (including photutils and scipy). The numerical
helpers use astropy FITS/WCS, numpy and photutils; plots use matplotlib. Download
requires public IRSA access. No credentials or account are needed.

```powershell
$skySource = 'tools/sky-plates'
$skyPython = 'python'
$skySpx = 'results/sky-plates-j0754-spherex-20261108'
& $skyPython -B "$skySource/fetch_spherex.py" --help
& $skyPython -B "$skySource/extract_spherex.py" --science "$skySpx/science" --calibration "$skySpx/calibration" --out "$skySpx/local"
& $skyPython -B "$skySource/compare_spherex.py" --root $skySpx --models 'results/sky-plates-j0754-models-20261108/models' --dust 'results/sky-plates-j0754-dust-20261108'
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

`fetch_spherex.py` preserves bounded native cutouts and per-request provenance,
reuses cached files, and retains prior failed receipts when retrying. Official
calibration files and their UUID/MD5 manifest must be obtained separately: see
the result directory's `method_audit.txt` and `calibration/` receipts. Calibration
selection in the extractor is explicit by detector and spectral UUID. Keep
alternative calibration versions in nested directories to avoid ambiguity.

The extractor samples the current pixel-convolved ePSF on native pixel centers,
normalizes on an enclosing stamp before cropping, and fits an inverse-variance
point source plus local plane at the fixed catalogue position. SOURCE pixels
are retained for the target and excluded from background fitting. Primary SIP
distortion is isolated from the alternate linear parent-pixel WCS. Wavelength,
bandwidth and solid angle come from matched native calibration maps. No
brightness or S/N clipping is applied. Flags, failures, negative measurements,
quality decisions and source hashes are retained. Aperture, mask, background and
coarse/fine centroid variants are diagnostics, not alternate primary measurements.

The official public IRSA Spectrophotometry Tool was also run for the same
coordinates, yielding 223 native rows in `native_IRSA_photometry.tbl`. This
server output must be saved separately before running `compare_spherex.py`.
The local extraction retains 220 fits, 182 under the declared primary quality
rule, and three failures. Local and IRSA reductions share science pixels and
are not independent observations. IRSA's D5/D6 wavelengths use the newer solution
(QR2 CWAVE +0.0064 microns), verified against the QR3 map; the local spectrum
preserves the science header's original QR2 solution. Native flux errors contain
a 3% quadrature term in addition to Tractor errors, verified numerically.

Several-percent reduction differences, high short-wavelength fit residuals and
fine-centroid sensitivity prevent precision seasonal-variability claims from
the formal local covariance errors. The nine-point ±0.3-pixel centroid grid is
too coarse to establish subpixel stability; the final extractor also runs a
bounded ±0.15-pixel continuous centroid sensitivity fit and retains success and
boundary indicators. The independent review's finer grid is retained with the
results. A generic detector-edge normalization bug
was caught with a failing regression test and fixed; the present accepted
cutouts happened to include all nonzero ePSF support, so it did not alter their
fluxes materially.

The SED check fits only 0.75–2.4 microns, then compares longer wavelengths. It
uses approximate rectangular CWAVE/CBAND responses, coarse CK04 atmospheres,
static Draine foreground attenuation and free normalization. Scores, best
parameters and a post-hoc 5% log-RMS sensitivity set are exploratory grid
results, not likelihoods, parameter measurements or confidence intervals.
An asynchronous/static spectrum cannot reconstruct the missing 2022 spectrum
or validate a period. No dominant warm infrared emitter is required by this
check; weak excesses, cold material and small line-of-sight clouds remain open.

Archive and literature searches recovered GALEX NUV, Swift UVOT and a plausible
Swift XRT counterpart, saturated UHS images, sparse mostly saturated PTF data,
and a 2011 direct-target conference-poster atmosphere fit. Association, optical
red leak, saturation, missing limits and literature status are recorded in the
result reports. They are not new-discovery or mechanism claims.
