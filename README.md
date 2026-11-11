# J0754 stellar variability laboratory

**Tang's other two stars, 11 November:** [full analysis](results/tang-sibling-analysis-20261111/report.txt)
and [main figure](results/tang-sibling-analysis-20261111/downlink_offset_models.png).
J0830 has a strong 38–39-day TESS modulation, close to its 38.389-day Gaia
candidate period, after allowing separate QLP reference offsets at downlinks.
J0736 favors roughly 19–21 days. J0830's independent 2015 K2 epoch does not
repeat a dominant 38-day signal. No excursion passes the initial seven-file
screen. All seven previously downloaded products are now analyzed; one K2SC
alternative is also retained. Candidate orbits and the historical fading cause
remain open.

**Comparison-giant screening, 11 November:** [report](results/giant-screening-20261111/report.txt),
[full catalogue screen](results/giant-screening-20261111/screened_catalog.csv),
[five physical matches](results/giant-screening-20261111/physical_shortlist.csv) and
[Kepler curves](results/giant-screening-20261111/kepler_screening_lightcurves.png).
The 61 published giant-labelled entries were screened; twelve Kepler light curves
were reduced for five selected stars. KIC 2852961 is the strongest activity
benchmark. Tang's other two historical fading stars have Gaia candidate SB1
periods near 38–39 days. That screening snapshot retains seven sibling TESS/K2
products; their completed follow-up is linked above. The Dixon 7,286-star rows
were not retrieved or screened.

**Current evidence checkpoint, 11 November:** [source-linked evidence table](docs/J0754_evidence_2026-11-11.csv)
([readable HTML](docs/J0754_evidence_2026-11-11.html), [text summary](docs/J0754_evidence_summary_2026-11-11.txt))
consolidates 36 entries across 37 retained result families. The
[common four-event TESS report](results/sky-plates-j0754-flare-comparison-20261111/report.txt)
adds durations, conditional energies, matched image injections and recovery
tests. All four excesses stay positive; baseline uncertainty and unresolved
hourly position diagnostics are retained explicitly.

**11 November follow-up:** the official eROSITA Gaia/NWAY catalogue selects
J0754 with p_any=0.919. Corrected GALEX fluctuations remain a weak, bin-dependent
UV variability hint. See the [report](results/sky-plates-j0754-galex-erosita-20261111/report.txt)
and [calibrated UV plot](results/sky-plates-j0754-galex-erosita-20261111/galex_corrected.png).

Reproducible archival astronomy tools and supporting evidence for **DASCH
J075445.9+164141**, an active giant with historical and modern fading episodes.
This repository unifies the plate photometry, modern optical/infrared analysis,
TESS flare and rotation work, spectroscopy, physical models and unsuccessful
tests into one research codebase.

**Earlier narrative: [10 November findings report](docs/J0754_findings_2026-11-10.txt)**
([illustrated Word version](docs/J0754_findings_2026-11-10.docx)). It distinguishes
measured results, model-dependent quantities and open hypotheses. Read the newer
evidence checkpoint above for subsequent measurements and limitations. This work has
not been peer reviewed. The star and its historical variability were reported
by [Tang et al. (2010)](https://arxiv.org/abs/1001.1395).

**10 November follow-up:** an expanded TESS search finds three additional
source-centered flare-like candidates (16 and 22 November 2021; 28 November
2023), alongside the original December 2021 event. Five new calibrated pixel
cutouts, the complete search hits, injection trials and rejected excursions
are included. See the [search report](results/sky-plates-j0754-flare-search-20261110/report.txt)
and [four-event comparison](results/sky-plates-j0754-flare-search-20261110/four_event_lightcurves.png).

**Multiwavelength inventory:** all seven follow-up archive groups were queried.
GALEX individual photons, two Swift UVW1 image extensions and three VLASS radio
cutouts were recovered. The official eROSITA DR2 Gaia counterpart catalogue selects J0754
(see the 11 November follow-up above). Current XMM/Chandra source searches and APOGEE DR17/DR19 target
queries are empty. See the [instrument inventory](results/sky-plates-j0754-multiwavelength-20261110/report.txt).

**Radio/UV measurements:** no radio detection or significant epoch change is
found in the three 3-GHz images. Approximate image-level 3-sigma bounds are
0.43/0.54/0.38 mJy (2019/2021/2024). UVOT exposure maps explain the first image's
low target counts; the two GALEX visit means are consistent. See the
[measurement report](results/sky-plates-j0754-radio-uv-20261110/report.txt) and
[radio comparison](results/sky-plates-j0754-radio-uv-20261110/radio_three_epochs.png).

## What is here

| Area | Tools and methods |
| --- | --- |
| Photographic plates | Comparison-star selection, quality cuts, same-exposure differential photometry, pixel/background fitting, scan-direction checks and historical coverage tests |
| Optical and infrared | ZTF, ASAS-SN, legacy ASAS and NEOWISE reductions; camera/filter separation; matched-state colors and control tests |
| Rotation | Stable-waveform checks, independent seasonal phases, period/alias searches, amplitude-brightness comparisons and shorter-window fits |
| Flares and eclipses | Calibrated TESS pixel localization, flux integration and passband energy corrections, normalized optical profile comparisons, multiday eclipse searches and injections |
| Spectra | FAST line/EW covariance propagation, registration, blueshift and RV repeatability checks; Gaia XP and SPHEREx shape comparisons |
| Physical hypotheses | Stellar-atmosphere and two-temperature spot models, extinction and partial covering, composite spectra, cloud budgets, recurrence and seismic diagnostics |

The main scientific findings are a replicated approximately 0.31-mag modern
optical fade and recovery, smaller infrared fading, a localized 2021 December
23 TESS superflare candidate, and recurring roughly 36-day modulation with
changing seasonal amplitude and phase. A Gaia candidate SB1 period is close to
that timescale. Evolving surface activity is the leading working explanation;
orbital confirmation and the cause of the long fade remain open.

## Install and test

Use Python 3.14, matching the original numerical environment. From this
repository root:

```bash
python -m venv .venv
# Activate .venv using your shell's usual activation command.
python -m pip install -r requirements.txt
python -B -m unittest discover -s tools/sky-plates -p 'test_*.py'
```

The current snapshot passes **205 tests** (147 in the initial release). Tests cover model
integration, injections, nuisance projection, period/alias behavior, covariance,
pixel reductions and failure cases; they do not substitute for independent
scientific review. No credentials are required for these offline tests.

## Reproduce a result

Run scripts from the repository root so their `results/` paths resolve:

```bash
python -B tools/sky-plates/search_seasonal_rotation.py
python -B tools/sky-plates/compare_rotation_brightness.py
python -B tools/sky-plates/check_short_rotation.py
python -B tools/sky-plates/check_rotation_aliases.py
```

These use the included daily/differential photometry. They regenerate outputs
in the corresponding result directories, so use a branch or copy if preserving
the publication snapshot matters. Larger image reductions require the original
archive products described below. Per-component methods, run commands and
limitations are in the [research guide](tools/sky-plates/README.md) and its linked
specialist READMEs. Report generation additionally needs
`requirements-report.txt`.

## Supporting evidence and original observations

`results/` retains scientific CSV tables, structured outputs, protocols, reports,
figures and selected compact public archive products, with the original stage
directory names to keep script inputs reproducible. Earlier exploratory and
failed reductions remain visible; the current report takes precedence over
superseded figure labels. A failed fixed waveform is not the same hypothesis
as rotation with changing seasonal phase.

The complete native FITS/calibration/image collection is too large to place in
ordinary Git history. [The data manifest](docs/native_data_manifest.json) records
paths, sizes, original SHA-256 hashes and which products are included.
[Public archive request URLs](docs/public_archive_requests.json) preserve stable
retrieval routes. The archive fetch scripts remain available; ASAS-SN job access
links, browser state and account receipts have been excluded. Retrieval depends
on archive availability and may require repeating a public query. This is not
a claim that every large reduction can run without downloading its inputs.

## Publication and provenance

This standalone release was exported from the astronomy component of a larger
private workspace on 10 November 2026. Only the astronomy snapshot is published;
the original workspace and its Git history are not part of this repository.
Personal account paths and private contact addresses were removed, document
authoring metadata was cleared, and credentials/private operational receipts
were excluded. See [the publication audit](docs/PUBLICATION_AUDIT.md).

Prepared by Kyle Hicks with AI-assisted analysis and coding by Ari, with
contributions from Gemini Argon. Astronomical observations retain their original
survey attribution. Please cite the original surveys and relevant methods when
using their data, as well as this repository when using its software. Open an
issue to report a reproducibility problem or propose another physical test.

Original software is released under the MIT license. Third-party observations,
atmosphere grids and calibration products retain their providers' terms and
required acknowledgments; the software license does not relicense those data.
