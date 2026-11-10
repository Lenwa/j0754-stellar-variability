# J0754: tidal synchronization and activity plausibility

9 November 2026. Literature and conditional scale assessment, not an orbital confirmation or stellar-evolution fit.

## Main finding

The proposed link is established astrophysics, not an improvised explanation: close giant binaries can be spun up and magnetically active, and tidal configuration may enhance activity even relative to single giants with similar rotation. J0754's candidate almost circular 35.92-day orbit and independently recovered roughly 36-day photometric modulation therefore make a coherent tidal-activity hypothesis. The available information does not measure the tidal torque, establish the prior spin history, or show that tides specifically caused the long dimming episode.

## What the literature adds beyond the existing local reports

The earlier `mechanism_assessment.txt` was read before this search. It already discusses Olah et al. 2013/2014 active K-giant analogs and their explicit comparison to Tang's objects. Those remain useful analogs, not new J0754 observations.

[Gehan, Gaulme & Yu 2022, A&A 668 A116](https://arxiv.org/abs/2211.01026) compared chromospheric indicators for thousands of Kepler giants. Their result is unusually pertinent: close-binary giants in spin-orbit resonance have stronger chromospheric activity than single giants of comparable rotation and physical properties. Thus the companion could matter through more than simply supplying a fast spin. This supports the physical mechanism in a population; it does not classify J0754 or prove a magnetic origin for its dimming.

[Gaulme et al. 2020, A&A 639 A63](https://arxiv.org/abs/2004.13792) found rotational modulation in 370 of approximately 4500 Kepler giants. Of 85 followed spectroscopically, 34 were clear binaries; 26 of 30 active nonoscillators were binary. In contrast, about 85% of active giants with detectable oscillations were not close binaries. The split matters: these figures must not be converted into a companion probability for J0754, whose analogous seismic classification is not established. Active single giants can reflect mass/evolution-dependent retained angular momentum, or previous engulfment/merger events. An active star does not require a presently detected companion.

[Dewberry & Wu 2025, ApJ 984 137](https://arxiv.org/html/2501.13929v2) compare tidal models with roughly 30,000 Gaia giant binaries. Theory predicts circularization around separation/radius 3–4, while some less evolved giants show circularization extending farther and circular systems out to about 10 radii. The model/data tension cautions against treating a single analytical timescale as definitive. At the illustrative J0754 primary mass 1.5 solar masses and minimum companion mass 0.2404, a/R is 6.06 for R=9.087 and 4.53 for R=12.16 solar radii: broadly relevant to this regime, not proof of synchronization. Their sample requires Gaia significance >15; J0754's 13.31 would fail this cut. Population context is not an independent validation of this candidate.

## Verified formula and conditional calculations

[Hurley, Tout & Pols 2002, MNRAS 329 897](https://arxiv.org/html/astro-ph/0201220v1), equations 25–33, supplies the conventional convective equilibrium-tide prescription. Native HTML, PDF and retrieval/hash receipt are retained. Using q=M2/M1 and beta=I/(M1 R²), the local near-circular synchronization inverse timescale is 3(k/T)q² beta^-1 (R/a)^6. The synchronized circularization inverse timescale is (21/2)(k/T)q(1+q)(R/a)^8. k/T=(2/21) fconv (Menv/M1)/tauconv; tauconv=0.4311[Menv Renv(R-Renv/2)/(3L)]^(1/3) years in solar units. The reduction factor is min[1,(Ptid/2tauconv)²], where inverse Ptid is abs(inverse Porb minus inverse Pspin). These are model prescriptions, not measured dissipation coefficients.

All following quantities are conditional. Use candidate P=35.9204566 days and f(M)=0.00458720 solar masses; an assumed primary mass 1.5 solar masses and edge-on orbit give M2=0.240409 and a=55.1088 solar radii. Adopt L=54.64 solar luminosities, envelope mass 0.7 M1, envelope depth R and beta=0.1. These structural values are explicit examples, not inferred stellar properties.

- R=9.087: convective turnover 101.07 days; near-synchronous local spin e-folding approximately 0.268 Myr; circularization approximately 3.9 Myr.
- R=12.16: turnover 122.73 days; spin e-folding approximately 0.0567 Myr; circularization approximately 0.46 Myr.

Both selected photometric period peaks, 35.78 and 36.30 days, imply long beat/pumping times and fconv=1 for these examples. A star initially rotating much more slowly need not evolve on these same short times: fconv becomes smaller far from synchronization. For the same example structures, changing the instantaneous spin to 100 or 300 days lengthens the local spin timescale, without yielding an integrated time to synchronize. For the R=9.087 example, the local spin times are 3.48 Myr at a 100-day spin and 6.58 Myr at a 300-day spin; for R=12.16 they are 1.09 and 2.05 Myr. These remain instantaneous values, not integrated synchronization durations. See tidal_scalings.json for exact values and equations.

The numerical grid varies assumed M1=1/1.5/2.6, R=9.087/12.16, inclination=90/60/30 degrees, envelope fraction=.3/.7/.9, envelope depth=.5/1 times R and beta=.05/.1/.15. Near-sync spin times span roughly 0.00185–3.62 Myr across this deliberately broad assumed grid. This is sensitivity, not a credible interval, an exhaustive structural range, or proof that synchronization has already occurred. These choices are not calibrated to J0754. No arbitrary tidal Q was inserted.

Synchronization and circularization are different. In this approximation the former scales with (a/R)^6 and the latter with (a/R)^8, with additional mass/inertia factors. A nearly circular orbit can be consistent with tides but can also reflect initial conditions or previous evolution; it does not measure a present synchronization timescale. The calculations hold stellar structure and orbital parameters locally fixed, omit magnetic braking/structural evolution and assume aligned spin/orbit. The star's residence time at this radius, envelope structure, inertia and evolutionary state are unknown. A forward stellar/binary evolution model is needed before a measured-age comparison is meaningful.

## Radius and interpretation boundary

The companion investigation's radius_consistency.json explicitly shows inconsistent atmospheric/evolutionary inputs. Fixed L=54.64 and T=4500–5206 K imply R about 12.16–9.09 solar radii; the poster logg=1.7 combined with assumed M=1–2.6 implies far larger radii, some close to their conditional Roche lobes. Do not merge these into one precise star. This report uses the luminosity/temperature branch only as a conditional calculation, not as an established radius. Near-Roche geometries would also invalidate a simplistic weak-distortion detached-star interpretation.

## What can actually be concluded

The tidal mechanism is physically plausible on ordinary giant evolutionary timescales under reasonable illustrative structures. The period match, activity and circular candidate orbit fit that picture, and newer giant-population evidence makes it more than a loose analogy. But the causal chain still has separate unmeasured links: independent orbit confirmation, real rotation/broadening and radius, tidal synchronization rather than coincident activity-related RV, then spot/activity changes causing the years-long flux decline. A star can be active without present tidal forcing; a close binary can be active without tides explaining a particular slow dip. No probability or unique cause follows from these timescales.

## Preservation

Saved: tidal_scalings.json, tidal_formula_source_receipt.json, Hurley2002.html/pdf, this tidal_literature_report.md. Formula and parameter values in JSON permit reproduction. No maintained source modified, no accounts, contact, payment or observing request. The PDF downloaded successfully; an initial optional text extraction failed because PyMuPDF was absent, so the native HTML supplied verified equations instead. No dependency installed. Parent owns combined report and Workspace continuation.


