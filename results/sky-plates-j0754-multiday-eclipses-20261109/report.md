# J0754: multi-day eclipse test using existing TESS data

9 November 2026. This extends the earlier 1/3/6-hour event tests, which could not constrain multi-day eclipses after their one-day highpass.

## What the data now establish

In the saved QLP SAP extracted flux, a persistent simple 1–6-day eclipse of 3% aperture depth at the candidate 35.92046-day period is substantially larger than the observed eclipse-like structure and is recovered at all 24 injected phases. This survives linear versus quadratic sector baselines and period changes of plus/minus 0.030613 days. Thus deep, repeatable eclipses of this particular kind are disfavored by the existing extraction. This is not a formal confidence exclusion or a universal companion limit.

At 1% depth, recovery depends strongly on phase, duration and baseline model. At 0.1–0.3%, unflattened QLP SAP does not provide reliable recovery under this deliberately conservative criterion. A shallow small-dwarf transit remains open. These are observed-aperture depths: dilution, unknown upstream extraction losses, grazing events and different eclipse shapes alter the astrophysical translation.

## Actual search

All eight existing native FITS were read using the retained QUALITY=0, finite positive flux/error rule. No brightness clipping, flare masking or new data retrieval. Each sector is normalized to its native median. Native data and native-sampled templates are averaged into occupied one-hour bins. Unlike previous work, there is no one-day highpass.

The four product groups are QLP SAP, QLP corrected (KSPSAP/DET), TESS-SPOC SAP and TESS-SPOC corrected (PDCSAP). These reductions share photons and are not independent confirmations. A template bank uses 288 uniform phases at period 35.9204566 days, durations 1/2/4/6 days, and trapezoids with 20% ingress and egress. Depth is jointly fit with per-sector offsets and linear slopes; quadratic baselines form a sensitivity. Equal-hour weighting avoids adopting misleading formal white-noise errors. Both data and templates are projected against the same nuisance columns.

The largest positive search-score QLP SAP trough has 6-day duration and fitted depth 0.733% with linear baselines or 0.480% with quadratic baselines. Its selected mid-eclipse phase is 10.4768 days from BTJD 2500 modulo the candidate period. The best corrected QLP trough instead has 1-day duration and only 0.021% depth, at unrelated phases. This discrepancy warns that corrected data can remove broad signals.

TESS-SPOC SAP's largest-score trough has 6-day duration at phase about 26.2–26.3 days, with 1.178% or 0.491% fitted depth for linear/quadratic baselines. Its corrected results also change with the baseline model. These are descriptive troughs, not eclipse detections.

## Is the largest pooled QLP trough actually repeated?

A post-hoc check fixes its phase and 6-day shape and fits each sector separately:

| QLP SAP sector | Linear depth | Quadratic depth |
|---|---:|---:|
| 44 | 1.445% | 1.284% |
| 45 | 0.163% | 0.074% |
| 72 | 0.347% | -0.119% |

Sectors 46 and 71 do not sample that template's eclipse interval. The pooled structure is therefore not a stable repeated eclipse of constant depth. At that phase the SPOC result is supplied only by sector 72, and differs materially from QLP. The three sectors are not independent detections of one eclipse. Full coefficients and template energies are saved in selected_trough_per_sector.json.

## Injection recovery

Depths 0.001, 0.003, 0.01, 0.03 and 0.1 were injected at 24 uniform phases for each duration. These enter the already extracted, normalized flux, not raw pixels or the upstream pipeline. Additive relative-flux injections approximate multiplicative aperture dimming; small native baseline variations are not multiplied by the injected signal.

After injection, the complete phase/duration bank is searched again. Success requires the winning center within half the injected duration, duration within a factor of two, recovered depth between 0.5 and 2 times the injected depth, and score above the largest positive score found in the original data. The last threshold measures whether the injection exceeds the actual existing structured variations; it is not a calibrated noise-only false-alarm threshold. Linear projection permits exact calculation of every injected fit without discarding the observed correlated residuals.

Recovered phase counts out of 24, in duration order 1/2/4/6 days:

| QLP SAP depth | Linear baseline | Quadratic baseline |
|---|---|---|
| 0.1% | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| 0.3% | 0 / 0 / 0 / 0 | 0 / 0 / 1 / 3 |
| 1% | 1 / 8 / 18 / 18 | 19 / 24 / 24 / 24 |
| 3% | 24 / 24 / 24 / 24 | 24 / 24 / 24 / 24 |
| 10% | 24 / 24 / 24 / 24 | 24 / 24 / 24 / 24 |

All 288 search phases have some QLP in-eclipse coverage: for a 1-day event there are 39–100 occupied hour bins aggregated across sectors, increasing for wider events. This does not establish complete individual eclipses or continuous phase completeness. SPOC has actual phase gaps for 1/2-day events; even 10% injections are not recovered at every tested phase. No absence claim should be based on SPOC alone.

Corrected QLP recovers even 0.1% artificial post-pipeline signals in this experiment, but that is NOT a physical 0.1% detection limit: a real broad signal could have been attenuated before these injections. The disagreement with SAP and SPOC is itself informative.

## Period sensitivity and scope

After the initial run, QLP SAP was repeated at 35.8898435 and 35.9510697 days, the candidate plus/minus one reported formal period error. At each period and both baseline orders, 3% injections recover at all 24 tested phases for all durations. The selected existing trough remains 6 days with depths about 0.45–0.77%, moving in fitted phase. This is a bounded sensitivity, not marginalization over an orbital posterior or a test of arbitrary periods.

The results constrain persistent trapezoidal aperture dimming at these trial periods, not intermittent occultation, every inclination, every companion radius, or the cause of the years-long decline. They neither confirm the Gaia orbit nor establish that the measured 36-day modulation is orbital. No more sensitive universal companion exclusion follows.

## Reproduction and verification

Maintained source: tools/sky-plates/multiday_eclipse_check.py and test_multiday_eclipse_check.py in build/github-desktop-review. Run from outer Lumen using runtime/sky-plates/venv/Scripts/python.exe -B and the script path. Inputs are the unchanged FITS in results/sky-plates-j0754-checks-20261108; native_inventory.json retains hashes.

Saved protocol, results.json, period_sensitivity.json, selected_trough_per_sector.json and six_day_profiles.png. Four new unit tests verify periodic trapezoid shape, joint polynomial-nuisance recovery, linear injection increments and blind recovery at a nondefault period. The full component suite passed 134 tests at this checkpoint. Source checkpoint accompanies the final handoff. Figure inspected. No new observations, downloads, contacts, paid services or publication.

