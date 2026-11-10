# J0754 stellar variability laboratory

Reproducible archival astronomy tools and supporting evidence for **DASCH
J075445.9+164141**, an active giant with historical and modern fading episodes.
This repository unifies the plate photometry, modern optical/infrared analysis,
TESS flare and rotation work, spectroscopy, physical models and unsuccessful
tests into one research codebase.

**Start with the [findings report](docs/J0754_findings_2026-11-10.txt)**
([illustrated Word version](docs/J0754_findings_2026-11-10.docx)). It distinguishes
measured results, model-dependent quantities and open hypotheses. This work has
not been peer reviewed. The star and its historical variability were reported
by [Tang et al. (2010)](https://arxiv.org/abs/1001.1395).

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

The publication snapshot passes **147 numerical tests**. Tests cover model
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
