# J0754: century-scale dimming pattern

## Result
The archive supports repeated, unequal multiyear fading and recovery, rather than a demonstrated clock. The deep 1930s episode, a smaller mid-1940s event, the 2006–2007 trough and the modern minimum around 2022 are useful landmarks. A smaller 2016 trough is also present in one well-sampled ASAS-SN V camera. The apparent approximately 15-year spacing of selected landmarks is worth recording, but neither it nor an approximately 44.5-year recurrence is established.

## Data and method
This analysis uses only the saved, previously screened DASCH differential photometry, target and matched-control plate curves, legacy ASAS exact-frame differential annual summaries, and the saved ASAS-SN light curve. Source paths and SHA-256 hashes are in results.json. There are 3,288 catalogue exposure rows, 854 detections, 583 primary-quality accepted detections and 552 retained differential measurements. These are different denominators; exposure coverage is not a calibrated faint-source detection-completeness model. Rejected or missing detections were not assigned a brightness.

Annual points are median magnitudes at the median actual observation date, not presumed dates of the true minimum. Historical series and modern camera/filter streams are kept identifiable. Modern annual-minimum comparisons require at least ten distinct observing days. The modern plot uses a separate baseline for each stream; connecting lines guide the eye through annual summaries, not measurements between them. No cross-passband absolute magnitude alignment or combined century-long periodogram was used.

The protocol was recorded before these calculations. The additional approximately 15-year hypothesis was explicitly recorded as a retrospective amendment before testing it. All event choices and recurrence hypotheses remain exploratory: no false-alarm probability or period confidence interval is claimed.

## Historical chronology and controls
- The early 1892–1893 feature has only three retained detections. It cannot supply a reliable duration or cycle epoch.
- The major episode is already faint in late 1931, with the faintest annual median in 1933 (14.113 mag at median observation date 1933.531). Annual medians recover to approximately 13.26 mag in 1941–1943. The approximately 0.85 mag annual-summary change is descriptive in the pooled historical system; it is not a universal cross-passband event amplitude. Sparse points in 1928 and 1931 leave ingress poorly bracketed.
- A smaller mid-1940s depression is present: pooled annual magnitude is faintest in 1946, not an arbitrarily chosen 1947. Within the ac plate series, equal-year means in 1941–1944, 1945–1948 and 1949–1956 give a 0.244 mag event excess above the linearly interpolated bracketing baseline. Restricting every annual median to at least three detections gives 0.246 mag, with 3/4/2 usable years and 16/18/9 measurements in the three windows. Nine matched-control contrasts under the default annual rule range from -0.115 to +0.125 mag; the target is larger than all nine. This supports a target-specific smaller event, but is not a calibrated significance test. The rh result of approximately 0.192 mag relies on sparse annual samples and fails the stricter coverage requirement, so it is supporting context rather than robust independent confirmation.
- A broad late-1970s/1980s change is visible, including fainter annual values in 1984–1985. For dnb, comparing 1969–1975, 1982–1986 and 1987–1989 gives a 0.079 mag excess above a linear bracketing baseline (0.080 mag with at least three detections per year). Several matched controls have comparable or larger contrasts. This feature is not a secure independent stellar cycle. The interval should not be described as perfectly flat.

The longest retained-detection gap is 1956.175–1964.131, nearly 7.96 years. Additional gaps include 1922.889–1928.293 (5.40 years), 1893.151–1897.241 (4.09 years), and 1928.293–1931.927 (3.63 years). No brightness state is inferred through these gaps.

## Modern chronology
Legacy ASAS V differential annual medians rise from 12.153 in 2004 to 12.308 in 2006 and 12.300 in 2007, then recover to 12.186 in 2009: a broad 2006–2007 trough, not an exact single-year timestamp. ASAS-SN camera ba in V has a smaller trough around 2016: 12.1905 in 2012, 12.312 in 2016 and 12.1895 in 2018. This is a within-camera descriptive feature, not a new comparison-star-controlled discovery.

All five usable ASAS-SN g camera streams have their faintest well-sampled annual bin in 2022. For example, bq has annual medians 12.643 in 2018, 12.9465 in 2022 and 12.758 in 2026. Recovery is partial. This corroborates the modern event timing; the previously analyzed ZTF and NEOWISE measurements provide additional wavelength information rather than another century-scale recurrence test.

## A bounded approximately 15-year hypothesis
Dividing the separation between the median dates of the 1933 and 2022 annual bins by six constructs a 14.7997-year interval. Agreement with those two anchors is guaranteed, not independent evidence. Predicted intermediate centers are approximately 1948.3, 1963.1, 1977.9, 1992.7 and 2007.5.

The mid-1940s feature and 2006–2007 trough lie near two of these centers, with unequal depths and broad timing. Around 1963 there are only two retained detections in the plus/minus three-year window. Around 1993 there is only a single 1989 edge point, not meaningful coverage of the predicted event. In contrast, the plus/minus three-year window around 1977.9 contains 68 retained measurements and does not display a replay of the deep 1930s event. The 2016 trough is an additional feature outside this proposed schedule. The catalogue ends before modern epochs, so zero historical counts for 2007 and 2022 mean outside DASCH, not absence of modern observations.

A schedule with freely varying depths and arbitrary missing cycles could accommodate these landmarks, but would not be a demonstrated predictive period. The observations favor episodic, unequal behavior over a simple fixed-depth repeating waveform.

## Explicit historical-waveform replay
The 1932–1943 annual recovery curve was linearly interpolated only inside its measured domain, shifted to the late-century observations, and tested with separate plate-series offsets. Missing ingress was not invented. A shift of 44.5 years samples 127 actual measurements spanning approximately 1976.8–1987.3. This full recovery-template domain differs from the plus/minus three-year event-center coverage above.

With series offsets only, the fitted template multiplier is -0.205. The residual RMS is 0.158 mag for no repeat and 0.372 mag for a fixed equal-depth repeat. Adding a free linear drift changes the comparison materially: no-repeat RMS is 0.148 mag and fixed-repeat RMS is 0.159 mag. The freely fitted multiplier is then -0.614 and is strongly entangled with the drift. Neighboring shifts of 43 and 46 years likewise do not yield a positive best-fit replay in the offset-only test.

Thus a fixed equal-depth, fixed-shape 1930s replay is not observed. The strong offset-only penalty is not a general exclusion of physical recurrence: an allowed long-term drift absorbs much of the broad recovery. Plate-response differences, censoring and uncertain baselines also limit interpretation. The template tests are deterministic sensitivity comparisons, not orbital or magnetic-cycle period measurements.

## Reproduction and files
Run decadal_patterns.py from tools/sky-plates using the existing sky-plates Python environment. It reads the saved inputs and writes results.json and decadal_pattern.png to this results directory. The source has four focused unit tests in test_decadal_patterns.py, alongside the existing component suite. No downloaded data, altered source measurements or new observations were used. protocol.md preserves the hypotheses and amendments; results.json includes exact numerical tables, annual counts, coverage windows, matched-control results, fit ranks and condition numbers.
