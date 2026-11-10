# Historical sky plates

## J0754 review report

`build_review_report.py` produces the compact scientific-review packet from
`j0754_report_content.json` and the existing outer `results/sky-plates-j0754*`
analyses. Run it from repository root with the bundled document Python providing
ReportLab and Pillow. It writes four narrative pages plus a figure appendix,
editable text, a complete original-resolution analysis PNG gallery, selected
JSON/protocol/report files, extracted CSV tables, research source and a hashed
inventory of retained native FITS. Outputs are under
`results/sky-plates-j0754-review-report-20261108`; native FITS and calibrations
are not duplicated into the compact ZIP. Published third-party paper/poster
pages and private account/job receipts are excluded. The packet is a review
snapshot, not a standalone raw-data archive or a publication submission.
Render the PDF with Poppler and visually inspect all pages after changes.

Experimental research scripts for checking long-term stellar variability against
comparison stars and independent plate archives. The initial target is DASCH
J075445.9+164141, reported by Tang et al. (2010), not a new discovery by this tool.

Maintained source: `tools/sky-plates` in the existing Lumen repository
`this repository`. Original downloads, protocols, environments and
generated figures belong in the repository root, outside this repository.
The first reconnaissance is `results/sky-plates-recon-20261107`; the authorised
pilot is `results/sky-plates-j0754-20261107`.

The pilot prediction and weakening outcomes were saved in `protocol.json` before
the light-curve reduction. The intended steps are catalogue-only comparison-star
selection, quality reduction, same-exposure differential photometry, calendar-year
block sensitivity checks, and comparison with APPLAUSE natural magnitudes and
plate images. Do not discard an astronomical signal merely because a variable
star was an outlier in a plate's calibration fit.

Use an isolated Python environment with `requirements.txt`. Local environment:
`python` in the outer workspace. Run numerical
checks with `python -B -m unittest discover -s tools/sky-plates -p test_analysis.py`
from this repository. Scripts use explicit input/output directories; commands are
documented here as they are completed.

## Reproduce the saved pilot

The following commands run from the **outer** Lumen workspace. The recon directory
must contain the retained original target light curve, APPLAUSE target extract and
exposure dates. Each download command reuses existing files; a new output directory
is required for a deliberate fresh archive snapshot. The astronomy environment was
tested with Windows Python 3.14. Its complete dependency freeze is retained with
the pilot output; `requirements.txt` pins the direct dependencies.

```powershell
$skyPython = 'python'
$skySource = 'tools/sky-plates'
$skyOut = 'results/sky-plates-j0754-20261107'
$skyRecon = 'results/sky-plates-recon-20261107'
& $skyPython -B "$skySource/fetch_controls.py" controls --out $skyOut --recon $skyRecon
& $skyPython -B "$skySource/fetch_controls.py" applause --out $skyOut --recon $skyRecon
& $skyPython -B "$skySource/analyze.py" --out $skyOut --recon $skyRecon
& $skyPython -B "$skySource/applause_diagnostic.py" --out $skyOut --recon $skyRecon
& $skyPython -B "$skySource/inspect_images.py" --out $skyOut
```

`controls` selects the nearest twelve catalogue-only matches (B within 1 mag,
B−V within 0.3, at least 300 detections and 20 arcsec catalogue isolation), and
retrieves the exposure/plate-series metadata. It starts at a 20 arcmin radius,
with the pre-specified 30 arcmin fallback if fewer than eight candidates exist.
Selection never uses the target's historical outcome or prunes controls by their
event contrast. The comparison ensemble is a robust same-exposure median, with
at least five controls per accepted target measurement.

The primary contrast is the median of annual medians in 1930–1939 minus that in
1910–1929 and 1950–1969; years require three observations. This is a descriptive
pilot statistic. The independent-year bootstrap does not include every correlated
systematic. Individual plate-series primary contrasts can be unidentifiable even
when a pooled or leave-one-series-out contrast is measurable. The separately
labelled sparse-series diagnostic permits one observation per year and must not
be described as primary within-series replication.

APPLAUSE source-table TAP queries timed out. The working fetcher uses the provider's
documented ten-arcminute cone-search API and filters desired plates offline.
`phot_range_flags=2` means a source is fainter than the estimated plate limit;
`phot_calib_flags=2` means calibration outlier, a different condition. All five
J0754 event rows in the first snapshot have the depth flag. Diagnostic handling of
those rows never overrides the pre-specified primary quality gate. Gaia G can be
missing in the plate table; match supplied catalogue synthetic natural magnitudes
and BP−RP colours in the same image instead. BP/RP calibration outputs from one
image do not constitute simultaneous independent colour measurements.

The plot leaves long observing gaps open. Scan displays are independently scaled
and are for checking gross defects, not for measuring a calibrated flux ratio.

The output is a bounded pilot. Plate defects, blending, colour response, source
splitting, poor limits and non-uniform sampling can change a result. A catalogue
binary-orbit fit does not establish the cause of historical dimming. Synthetic
tests validate calculations, not the archive measurements or astrophysics.

Primary documentation:
- https://dasch.cfa.harvard.edu/dr7/reduce-lightcurve/
- https://dasch.cfa.harvard.edu/dr7/colorterms/
- https://www.plate-archive.org/metadata/applause_dr4/source_calib/
- https://arxiv.org/abs/1001.1395

## Independent image follow-up (8 November 2026)

The next authorised experiment is retained in outer
`results/sky-plates-j0754-followup-20261108`. Its original protocol and measurement
specification precede pixel reduction. A dated design amendment records the move
from the incomplete extracted-source catalogue to complete Gaia catalogue-forced
measurements. The later isotonic calibration is explicitly secondary, designed
after the primary quadratic calibration failed. Both verdicts are **inconclusive**
under the saved screening rules. Exploratory estimates point in the reported
direction, but do not independently confirm the event or its physical cause.

Native FITS crops cover two 1919 and three 1935–1937 Hamburg Lippert Astrograph
plates. `fetch_applause_pixels.py` requires bounded HTTP byte ranges and refuses a
full-file response. It reads a complete row band (at most 90 MB per scan) and saves
only a 1300-pixel crop. The original headers, transfer receipts and checksums are
retained. No resampling of the image or injection of synthetic incident flux is
used. Scanned photographic density is calibrated empirically against colour- and
magnitude-matched Gaia DR3 stars, with proper-motion propagation and registration
against independent plate centroids. Missing or negative control signals are
retained in the output, with explicit reasons for unscored predictions.

From the outer workspace, reusing the variables in the earlier example:

```powershell
$skyFollow = 'results/sky-plates-j0754-followup-20261108'
& $skyPython -B "$skySource/fetch_applause_catalogue.py" --old $skyOut --out $skyFollow
& $skyPython -B "$skySource/fetch_applause_pixels.py" --old $skyOut --out $skyFollow
& $skyPython -B "$skySource/remeasure_applause.py" --out $skyFollow --recon $skyRecon
& $skyPython -B "$skySource/remeasure_applause.py" --out $skyFollow --recon $skyRecon --estimator isotonic
```

The first command rebuilds the union from retained VOTables and downloads only
missing catalogue snapshots; the other fetcher also reuses existing files. The
data selection is specific to this target, not a general plate-survey pipeline.
To run the numerical checks from the source repository:
`python -B -m unittest discover -s tools/sky-plates -p 'test*.py'`.

Read `pixel_results.json` as the primary image experiment and
`isotonic_pixel_results.json` as the secondary calibration. The star-stamp figures
use the same density stretch within each row; display scales vary between plates.
Held-out RMS describes only scored controls, not complete recovery. A near-normal
comparison star is not recovered at its expected position on the 1937 scan; whether
this is variability, a local plate effect or another problem remains unresolved.
Blank-aperture tail counts are descriptive, not a Gaussian false-alarm significance.
`blank_noise_only_bright_bound_mag` excludes calibration and other systematic
uncertainty and is not a validated stellar brightness confidence limit. No model
extrapolation is used, and failed or unavailable predictions are not nondetections.

APPLAUSE acknowledgements: the digitized images are provided by Hamburger
Sternwarte (Universitaet Hamburg). Funding for APPLAUSE has been provided by DFG,
Leibniz Institute for Astrophysics Potsdam (AIP), Dr. Karl Remeis-Sternwarte Bamberg,
and Hamburger Sternwarte. The original FITS headers retain the full provider
acknowledgements and CC0 licence. Gaia DR3 is an ESA Gaia/DPAC data product.

## Enlarged fields, orthogonal scans and modern observations

Outer `results/sky-plates-j0754-checks-20261108` retains the next authorised stage.
Both scan orientations of all five physical plates were downloaded as 2200-pixel
crops. The 1937 X crop is clamped to the actual scan boundary, with its original
pixel origin retained; no outside-image pixels are manufactured. Transfer caps
were explicitly raised to 150 MB per row band for this stage. The larger Gaia
field has radius 0.65 degrees and came from the Heidelberg Gaia service after
ESA returned HTTP 503. The unmodified earlier APPLAUSE union supplies registration
references, while the complete enlarged Gaia field selects photometric controls.

Control counts increase to 25–28. The broader-field target estimates change
substantially and all complete corroboration verdicts remain inconclusive. The
1937 Y target has no catalogue calibration; its colour coefficient can be recovered
from other matched stars. The 1919 March Y scan has neither a target calibration
nor registration references; the paired-scan colour response is explicitly assumed
for diagnostics and its registration gate fails. Such estimates are not usable
independent photometric confirmations. Raw absence at the two troublesome 1937
control positions persists in both orientations. Their positions are near scan
boundaries, but whether this reflects coverage, defects or historical variability
is unresolved. X/Y fits also use different control samples.

The modern snapshot queries J0754 plus all sixteen original 1937 comparison stars
within three arcseconds. Fourteen sources have actual public ZTF rows, including
the target; three return valid empty tables. J0754 has 575 raw and 470 measurements
passing the stated flag/error screen, spanning March 2018–October 2025. Eight TESS
products cover 2021, 2023 and January 2026; some are different reductions of the
same epochs. The source and observation metadata are not counted as separate new
observations. All 17 numerical tests passed for this stage.

From the outer workspace:

```powershell
$skyCheck = 'results/sky-plates-j0754-checks-20261108'
& $skyPython -B "$skySource/fetch_applause_catalogue.py" --old $skyOut --out $skyCheck --gaia-radius 0.65 --gaia-service heidelberg --registration-from $skyFollow
foreach ($skyOrientation in @('x','y')) {
    & $skyPython -B "$skySource/fetch_applause_pixels.py" --old $skyOut --out $skyCheck --size 2200 --maximum-band-mb 150 --edge-clamp --orientation $skyOrientation
    foreach ($skyEstimator in @('quadratic','isotonic')) {
        & $skyPython -B "$skySource/remeasure_applause.py" --out $skyCheck --recon $skyRecon --estimator $skyEstimator --orientation $skyOrientation
    }
}
& $skyPython -B "$skySource/fetch_modern.py" --previous $skyFollow --recon $skyRecon --out $skyCheck
& $skyPython -B "$skySource/summarize_modern.py" --out $skyCheck
& $skyPython -B "$skySource/modern_differential.py" --out $skyCheck
& $skyPython -B "$skySource/scan_comparison.py" --out $skyCheck
```

ZTF results keep object/filter zero points separate and retain large excursions.
An exploratory correction of the largest target data set requires at least three
other stars with exactly the same exposure/filter/field/CCD/quadrant, with at most
one control per star per frame. It retains 364 target measurements and an apparent
2018-to-2022 annual-median change of about 0.31 mag, followed by partial recovery.
This diagnostic was designed after seeing the modern plots. It does not rule out
source-specific detector, saturation, colour or reference-star effects and does
not connect the apparent trend to the 1930s event or a physical mechanism.
TESS curves are sector-normalised, retain quality-zero samples and leave long gaps
open; detrending and independent sector normalisation can suppress slow changes.
Modern stability does not establish stability in 1937.

## Independent ASAS-SN attempt, 8 November 2026

The completed Sky Patrol V1 recomputation corroborates the direction of the
modern ZTF fading and recovery. Three g-camera streams give 2018–2022 changes
of +0.313, +0.285 and +0.3035 mag. All five adequately sampled g streams brighten
by 2025. The earlier direct-V2 coverage failure below is an archive-access
outcome; the new V1 measurements resolve that gap for this test.

The pre-inspection protocol asks whether ASAS-SN reproduces the apparent ZTF
2018–2022 fading and subsequent recovery. Native public Sky Patrol V2 catalogue
and photometry Parquet blocks were downloaded using the
[official client protocol](https://github.com/asas-sn/skypatrol). Gaia IDs remain
integers or decimal strings throughout. Camera names are case sensitive: `ba`
is V, while `bA` is g. Filters and cameras retain separate summaries.

All seventeen requested stars have catalogue matches. The native photometry
block has 18,391 rows for fourteen comparison stars and no target measurements.
A single-target lookup and a 36-arcsecond master-catalogue cone both identify
the same target and return zero photometric rows. These empty responses provide
no test of variability. The formal verdict is insufficient target coverage.

The older [ASAS-SN variable-star page](https://asas-sn.osu.edu/variables/05a882b3-3477-5b94-b81b-17a186f6b4c1)
does supply 222 V-band observations from camera `ba`, January 2013–November 2018,
through its JSON export and embedded HTML. The CSV export is header-only in this
snapshot, so it is preserved but cannot supply the measurements. These published
rows do not include quality flags; their summary screens finite values, positive
flux/error and `0 < mag_err < 99`, without brightness clipping. Annual medians
require ten observations. They end before the 2022 and 2025 comparison epochs.

The catalogue labels J0754 `ROT:` (uncertain rotational variable) and lists a
912.3158778-day period. These are catalogue claims, not a period independently
validated here or an explanation for the modern or 1930s fading. The catalogue
reference is ASAS-SN variable stars Paper II, published in 2019. Full-history
coordinate photometry through Sky Patrol V1 is the next route; its prepared
6000-day aperture-photometry request required a human CAPTCHA before submission.
Kyle completed that check; Ari inspected the offered older precomputed curve,
then submitted the recomputation. Job `21eb35a3-8487-4560-a40d-be000fae5f7d`
was observed queued and then started. Its public status endpoint is the
computation URL plus `/status`; the site's own client normally polls it every
three seconds. This stage uses occasional reads while retaining the receipts.
The job completed at 21:04 UTC. Its
[result page](https://asas-sn.osu.edu/sky-patrol/coordinate/d2c2ad05-c9a6-4362-93ac-523c2584d094)
supplies 5,249 native rows (965 V, 4,284 g), February 2012–8 November 2026;
5,180 pass the finite/positive/error/limit screen. V1 supplies no explicit image
quality flag. Sentinel values and limits remain in the preserved CSV but are
excluded from detection summaries. No brightness clipping is used.
The alternate precomputed photometry record also has the same 222 epoch/camera
keys, ending November 2018, with magnitude differences at most 0.002 mag from
the variable-star record. It supplies no additional observing epochs. Its JSON
API paginates by 100 rows, while the HTML table paginates by 20; follow the API's
returned `next` links when saving the full record.

```powershell
$skyIndependent = 'results/sky-plates-j0754-independent-20261108'
& $skyPython -B "$skySource/fetch_asassn.py" --previous $skyFollow --out $skyIndependent
& $skyPython -B "$skySource/summarize_asassn.py" --out $skyIndependent --ztf $skyCheck
& $skyPython -B "$skySource/fetch_asassn_job.py" --job-id '21eb35a3-8487-4560-a40d-be000fae5f7d' --out $skyIndependent --wait-minutes 12
& $skyPython -B "$skySource/asassn_replication.py" --out $skyIndependent --ztf $skyCheck
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

Existing native files are reused. Use a new output directory for a fresh archive
snapshot or to retry an interrupted catalogue/block query. The figure uses
separate axes for published V magnitudes and the earlier ZTF g differential
result, with no absolute zero-point alignment. Five additional tests guard
missing-target coverage, absent comparison epochs, camera/filter case and
malformed-response cache handling and polling limits. Four further replication
tests guard limits/sentinels, camera separation, day grouping and sparse dates.
That stage's numerical suite had 26 tests; the period test below brings the total to 34.
If a refreshed query supplies target rows, the missing-coverage summarizer fails
explicitly so the planned per-camera replication analysis must be implemented.
The job reader performs GET requests only, retains timestamped status responses,
and saves the completed result page and its actual CSV export. It never
submits another computation. A pending exit can be resumed with the same job ID.
The waiting limit must be finite and between zero and sixty minutes; malformed
limits are rejected before any request, preventing accidental rapid polling.

The replication uses annual all-image medians within each camera/filter, with
at least ten observations on ten UTC dates in each compared year. Its 95%
sampling intervals resample complete UTC-date groups (3,000 draws), retaining
the grouping of consecutive exposures. All three adequately sampled fading
streams have intervals above zero. Five recovery streams have intervals below
zero. Error/seeing screens and common-month/season comparisons retain the same
direction where adequately sampled; insufficient subsets are reported. Detailed
CSV-specific choices followed preliminary medians and are explicitly recorded
as such, while the directional prediction preceded the data.

All fourteen available original comparison stars are retained. Their median
2018–2022 changes in the three principal cameras are approximately -0.0017,
+0.0007 and +0.0001 mag. These use V2 coadded/differently reduced data, so they
are annual diagnostics rather than identical-frame V1 corrections. V2 coverage
in 2025 is insufficient for the equivalent recovery control test. Bootstrap
intervals exclude calibration systematics and inter-day correlation; camera
streams can share systematics. The observational pattern is corroborated across
surveys, while cause, period and any connection to the 1930s event remain open.
The plotted 2026 point is explicitly a partial year through 8 November.

## Testing the catalogue's 912-day period

The longer record does not support a stable 912.3158778-day clock as the main
source of J0754's modern variability. Tests use g data from 2019 onward, separate from
the older published V fitting sample, with single constant camera offsets.
Daily medians are combined into fixed 30-day bins without brightness clipping,
year/cycle centering, phase adjustment by cycle or detrending before the test.
The primary samples are 299 ASAS-SN bins and 55 ZTF differential bins.

A single-harmonic clock explains only 3.4%/3.7% of camera-centered variance in
ASAS-SN/ZTF. Allowing two harmonics raises this to only 7.9%/8.1%. The flexible
clock increases 90-day blocked prediction error by 12.7%/15.3%, and increases
retrospective forward error by 18.6%/25.3% when trained before 2023. Either
harmonic model could meet the saved support criteria; adding flexibility cannot
invalidate a successful simpler model. Neither meets those criteria here.
Adjusting the period within ±10%, selected only on training observations, also
fails to give useful later predictions.

Two cycles are nearly five years, providing almost the same observing season.
One-to-one same-camera pairs within 15 days of that separation differ by a
median absolute 0.1485/0.1745 mag (75/14 pairs). Tightening to five days retains
median discrepancies of 0.146/0.158 mag. Correlated pairs are not treated as
independent significance measurements. At the actual sampling times, stable
0.30-mag simulated sinusoids are readily recovered under the declared 0.05-mag
noise model; this is conditional power, not a false-alarm probability.

A five-knot slow reference describes the data much better with the same number
of parameters as the two-harmonic clock. It is an interpolation reference and
has no identified forecast for future knots. Broad exploratory scans favor
longer scales but cover less than one cycle of their best-fitting sinusoids;
they do not establish replacement periods. A weak or evolving 912-day component
is not excluded. The source catalogue's classifier confidence is not a period
confidence, and its original search was limited to 0.05–1000 days
([Paper II, section 3.2](https://arxiv.org/html/1809.07329v2#S3.SS2)).

The test protocol predates fitting in this stage, while the broad modern trend
was already known. Forecasts are retrospective temporal tests, not blind new
observations. The legacy 222-point V sample is summarized separately and never
counted as independent period confirmation. Bin-width and quality checks agree.

```powershell
$skyPeriod = 'results/sky-plates-j0754-period-20261108'
& $skyPython -B "$skySource/analyze_j0754_period.py" --independent $skyIndependent --ztf $skyCheck --out $skyPeriod
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

`period_models.py` contains camera-offset fits, temporal validation and matched
phase comparisons. Eight new synthetic tests guard stable-signal recovery,
slow-drift rejection, passband separation, duplicate dates, unique phase pairs,
held-out camera identifiability, forecast leakage and model-selection behavior.
All raw inputs remain in earlier snapshots; the new directory preserves input
hashes, the protocol, retained bins, complete results, a figure and report.

The ZTF data are provided through the NASA/IPAC Infrared Science Archive; TESS
products are provided through MAST, with QLP and TESS-SPOC reductions identified in
their original FITS headers. Gaia source positions and IDs, requests, and raw
photometric samples are retained alongside the summaries and figures.

## Wavelength and mechanism comparison

J0754's modern fading and partial recovery also appear in native NEOWISE W1/W2
photometry. This is a test of competing explanations, with ordinary foreground
dust, photospheric changes and occulting/binary geometries remaining distinct
hypotheses. It does not establish that the historical and modern events have
the same cause or a fixed recurrence interval.

`fetch_dust_data.py` obtains public IRSA Gator snapshots for the target and the
original fixed sixteen Gaia comparison stars, plus target AllWISE/2MASS tables.
The native typed IPAC tables, exact requests and CSV conversions are preserved;
character quality flags retain leading zeros. Existing native files are reused.
TAP cone queries timed out in the original run; Gator supplied the same catalogue.

`analyze_j0754_dust.py` applies separately declared band-specific quality cuts,
chooses duplicate candidates by position before inspecting brightness, and
subtracts a per-source-centered comparison ensemble on the exact scan/frame.
At least three controls and five frames per visit are required. No target value
enters the ensemble or control placebos. Visits are separated by gaps exceeding
30 days; annual medians of visit medians and matched calendar-half-year changes
are retained. Colors use paired bands on the same exposure. ASAS-SN optical
matches retain camera identity and require three UTC dates within ±30 days;
the ±15-day check is mostly insufficient for matching both years and seasons.

Primary 2018–2022 changes are +0.11625 mag (W1) and +0.089375 mag (W2).
Both observing seasons, raw photometry, quality variants and leave-one-control-out
checks retain the fading. All sixteen controls remain eligible under identical
measurement cuts; fifteen have sufficient primary W1 contrast and sixteen W2.
Their placebo changes span -0.00875 to +0.003875 mag (W1) and -0.027 to +0.02925
mag (W2). Official NEOWISE calibration documentation describes seasonal offsets,
so small uncorrected W2 changes alone would be inadequate evidence.

For a 0.30-mag optical fade, three published Draine Milky Way grain models
predict only 0.013–0.015 mag in W1 and 0.007–0.009 mag in W2. This challenges
those simple foreground screens as the entire explanation. Fixed-area Planck
cooling examples predict larger infrared fading but also undershoot the
measurements; no causal model has been fitted. These are monochromatic
illustrations, not passband-integrated stellar atmospheres or radiative transfer.
Unusual grains, partial occultation, mixed mechanisms, starspots, temperature
and radius changes remain open. Static asynchronous 2MASS/AllWISE colors are
not a fading color curve; no warm excess by itself uniquely identifies a ring.

The declared conditional bootstrap resamples complete UTC-date groups within
the four actual 2018/2022 visits, keeping visits and ensemble corrections fixed.
Its intervals exclude calibration, comparison and inter-date correlation errors;
only two or three dates per visit are available. They are sampling diagnostics,
not physical-model confidence intervals. No calibrated significance is claimed.
NEOWISE ends in 2024 and cannot test the 2025–2026 infrared recovery.

`compare_j0754_events.py` displays the retained historical and modern records
without optimized shifts, time stretching or amplitude matching. The 1930s
episode is recorded; an early 1892–1893 feature has only three measurements.
The 1973–1982 record has 106 retained detections and must be confronted by any
shorter recurrence proposal. A spacing between episodes is not yet a period.
Independent historical APPLAUSE photometric corroboration remains inconclusive.

```powershell
$skyDust = 'results/sky-plates-j0754-dust-20261108'
& $skyPython -B "$skySource/fetch_dust_data.py" --previous $skyFollow --out $skyDust
& $skyPython -B "$skySource/analyze_j0754_dust.py" --previous $skyFollow --asassn "$skyIndependent/asassn_v1_lightcurve.csv" --out $skyDust
& $skyPython -B "$skySource/compare_j0754_events.py" --historical 'results/sky-plates-j0754-20261107/dasch_differential.csv' --asassn "$skyIndependent/asassn_v1_lightcurve.csv" --out $skyDust
& $skyPython -B -m unittest discover -s $skySource -p 'test_*.py'
```

Eight numerical tests guard band flags/errors, explicit cut sensitivity,
target-free exact-frame correction, frame identity, visit coverage, descending
extinction grids, Planck normalization and date-group resampling. The complete
capability now has 42 tests. Results also contain a literature/catalogue audit
and conditional calculations for Gemini's suggested orbital geometries.
Gaia's 35.92-day SB1 candidate falls below the processing paper's recommended
gold-sample efficiency/significance cuts; it is a testable lead, not a validated
cause. Source, raw archive data, generated results and Python runtime remain
separate under the existing project organization.

## Competing models and extended archive coverage

The [spectral-change comparison guide](stellar_models_README.md) documents the
passband-integrated stellar/spot, finite-depth partial covering and spectral
dilution models, additional-color forecasts and exact 35.92-day photometric
coherence test. The [historical extension guide](historical_extension_README.md)
documents newly recovered 2002–2009 ASAS photometry, unchanged early DASCH quality
audits and expanded APPLAUSE/Catalina/NSVS coverage checks. Current results,
spectral quality audits and SPHEREx image-coverage metadata remain separately
under outer `results/sky-plates-j0754-models-20261108/`.

The existing data admit several physically distinct spectral-change models.
Forecast color ranges are sampled parameter/tolerance sensitivities, not model
probabilities. Historical plate footprints and SPHEREx spectral-image overlaps
are not new calibrated target fluxes. The recorded historical/modern pattern
remains worth investigating without assigning it a fixed long period or cause.

The [SPHEREx guide](spherex_README.md) documents the newly extracted 2025/2026
infrared spectrum, native IRSA comparison, exploratory photospheric checks and
wider UV/X-ray/archive search. Calibration systematics and physical-model
ambiguities remain explicit; numerical checks accompany the extraction tools.

The [physical diagnostic guide](physical_tests_README.md) applies equivalent-width,
temperature-feature, TESS variability/event, CO-window and conditional binary/model
tests to preserved observations. Supplied theoretical examples remain distinct
from measured inputs; the tests do not establish a unique fading mechanism.
The complete component suite now contains 114 numerical tests.
