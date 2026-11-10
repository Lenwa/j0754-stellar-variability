# J0754: how strongly does the evidence support companion-driven activity?

9 November 2026. New checks of existing data and primary literature.

## Assessment

The companion-driven activity interpretation has become more compelling as a working explanation. Two new checks help: the approximately 36-day photometric rhythm is preferred over its simple yearly sampling aliases, and a conventional convective-tide calculation permits synchronization on relatively short evolutionary timescales under explicit example structures. Published red-giant research independently supports a connection between tidal resonance and enhanced magnetic activity.

This does not yet complete the causal chain. The strongest established portion is an active giant with a recurring approximately 36-day timescale. A near-synchronized binary is a coherent interpretation of that plus Gaia's candidate orbit. The claim that this companion causes the multiyear dip remains less directly tested. No probability has been assigned.

## 1. A new check on the photometric clock

Previous owner work already recovered seasonal peaks at 35.78 days in ASAS-SN and 36.30 days in ZTF, allowing changing phase and amplitude. We reused those exact data and nuisance models rather than claiming a new discovery.

The new bounded comparison tested Gaia's 35.9204566-day period against annual frequency aliases at 32.7041699 and 39.8383535 days, half/double periods, and the synodic lunar month. Equal-season fractional reductions in residual sum of squares were:

| Survey | 35.9205 days | 32.7042 days | 39.8384 days |
|---|---:|---:|---:|
| ASAS-SN | 0.309 | 0.150 | 0.106 |
| ZTF | 0.559 | 0.155 | 0.287 |

The candidate period beats both yearly aliases in 7 of 9 ASAS-SN seasons and 5 of 6 ZTF seasons. Removing any one season preserves a positive average advantage over the stronger alias in each remaining season.

Noiseless signals injected at the actual selected observing times recover all their signal at their true period, but only about 20–24% at its yearly aliases. The cadence can create secondary peaks, but it does not make those clocks equivalent. These are deterministic sensitivity checks, not false-alarm probabilities or proof against every instrumental process.

A 71.84-day rotation with two similar active longitudes could still have its strongest harmonic at 35.92 days. The different selected survey peaks are not a measured differential-rotation rate. Their offsets from Gaia's period are only approximately -0.40% and +1.06%, without calibrated uncertainties for the evolving spots.

## 2. Can a companion plausibly maintain rapid rotation?

Yes, conditionally. The saved Gaia SB1 solution has P = 35.92046 days, K = 10.7245 km/s and e = 0.023 +/- 0.056. Its mass function is 0.0045872 solar masses. With an illustrative primary mass of 1.5 solar masses and an edge-on orbit, the companion minimum mass is 0.2404 solar masses and separation is 55.11 solar radii.

Using the cautioned luminosity estimate of 54.64 solar luminosities with temperatures 4500–5206 K gives radii 12.16–9.09 solar radii. Synchronous equatorial speeds would be 17.13–12.80 km/s. These are predictions to compare with a real rotational-broadening measurement, not measured v sin i.

Using the convective prescription in [Hurley, Tout & Pols (2002), equations 25–33](https://arxiv.org/html/astro-ph/0201220v1), an example envelope mass of 0.7 times the primary mass, envelope depth equal to its radius and normalized moment of inertia 0.1 gives local near-synchronous spin e-folding times of about 0.057–0.268 million years over those radii. Corresponding circularization times are about 0.46–3.9 million years.

These are local fixed-structure examples, not the star's measured history or time taken to become synchronized. Rotation initially far from synchronization reduces effective convective dissipation. Unknown envelope structure, evolution, braking and inclination matter. A deliberately broad assumed grid spans approximately 0.00185–3.62 million years near synchronization; that span is not a confidence interval. Detailed formulas and assumptions are in tidal_literature_report.md and tidal_scalings.json.

The atmosphere estimates must not be combined indiscriminately. The poster gravity log g = 1.7 would imply a much larger radius at ordinary assumed giant masses, near the conditional Roche lobe. This is a warning that the available stellar parameters are not one consistent measured solution, not evidence that the star overflows its Roche lobe.

## 3. What the literature adds

[Gehan, Gaulme & Yu (2022)](https://arxiv.org/abs/2211.01026) find stronger chromospheric emission in red giants in tidally resonant close binaries than in single giants with comparable rotation and physical properties. This supports precisely the sort of physical connection proposed here, but does not identify J0754 as a member of that class.

[Gaulme et al. (2020)](https://arxiv.org/abs/2004.13792) also show why activity alone cannot settle binarity: active giants include single rapid rotators and stars potentially affected by previous engulfment or mergers. Their binary fractions depend strongly on whether oscillations are detected; J0754 has no equivalent validated classification. Those population percentages cannot be used as its individual binary probability.

[Dewberry & Wu (2025)](https://arxiv.org/html/2501.13929v2) compare giant-binary circularization with Gaia and find limitations in simple tidal prescriptions. J0754's illustrative separation/radius ratio of about 4.5–6.1 lies in a physically relevant regime, but the candidate fails that paper's significance selection. It is contextual evidence, not an independent confirmation.

## 4. The decisive unresolved link

Gaia gives a useful candidate orbit, but the source's significance 13.31 and efficiency 0.1889 fail the recommended gold-sample cuts. Its good fit and near-circular shape are supportive, not sufficient. Near-circular periastron degeneracy partly explains low efficiency; failing a selection cut does not prove a false orbit.

There are no released individual Gaia RV measurements or RVS line profiles for this source. As the [Gaia SB1 processing paper](https://arxiv.org/html/2410.14372v1) explains, fitting an orbital-shaped RV curve does not establish the physical origin of every RV variation. A rotation-period signal in line shapes can mimic an orbital signal. Matching photometric and RV periods therefore has two readings: synchronization if the RV is orbital, or activity affecting both.

The critical missing evidence is independent high-resolution Doppler motion across the candidate orbit, checked against line shapes and activity indicators. No such observation was requested here. The recovered FAST spectra and current broadband SED do not supply that test.

Even after confirming a synchronized binary, linking it to the long dip would still require a credible relationship between changing atmospheric/spot properties and the slow brightness evolution. Rotational modulation measures asymmetry, not total spot coverage. The existing activity, flare-like event, multiband fading and successful spot-model consistency check make that an important hypothesis, without proving every arrow in the chain.

## Reproduction and preservation

New cadence results and unchanged-data hashes: ../sky-plates-j0754-alias-check-20261109/report.txt and results.json. Maintained source: tools/sky-plates/check_rotation_aliases.py, test_rotation_aliases.py and rotation_aliases_README.md in the existing source repository.

Tidal inputs, numerical examples and formula source are retained in this directory. radius_consistency.json contains transparent radius/gravity/Roche calculations. spot_centroid_toy.json is an appendix-only mathematical missing-flux bound, not evidence confirming a binary: actual cross-correlation velocities and atmospheric line changes need not obey its restricted model.

No new observations, contacts, paid services, original-data changes or publication. The 121-test component suite passed, including three new alias tests. Independent review verified all 1,296 tidal-grid rows within 2.9e-9 relative and found no blocking mathematical issue. The alias source, tests and guide were committed as 8e099023d76b161aac8edaa0be3bdca1cd9abf96; unrelated edits were preserved and nothing was pushed. Research results remain local outside source Git.
