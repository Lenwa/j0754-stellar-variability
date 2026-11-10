# Broader dust test protocol

2026-11-09. New tests use only saved observations, atmospheres and responses. Read prior model_comparison.json and physical constraints first. No new observation, download, or fitting to a newly acquired datum.

Training anchors: March g=0.312, W1=0.1175 mag; retrospective check W2=0.08875. October g=0.3225,W1=0.115,W2=0.09 provides sampling sensitivity. W2 is known but excluded from optimizing parameters. A 0.02 mag per-training-band envelope is descriptive, not a confidence set.

Extend standard Draine Rv3.1/4/5.5 partial covering to nonzero initial coverage (0,0.1,0.3,0.5,0.7) with fixed optical depth and increasing coverage. Also vary optical depth at fixed coverage, with nonzero initial optical depth. Integrate every attenuation through the existing passbands and CK04 spectra. Compare nominal5000K and4500/5500K atmosphere sensitivity at [M/H]=-0.5,logg2.5. Do not infer physical grain sizes from a fitted extinction slope.

Fit phenomenological power-law extinction alpha=0–2.5 and gray-plus-standard attenuation, with g/W1 only. Explore partial power-law covering and dust-plus-spots, documenting nuisance parameter dimensions and predictions, not declaring additional freedom evidence. Mixed model uses fixed initial spots and final spot coverage plus a foreground screen; both components attenuated before band integration.

Retain sampled training-compatible ranges as well as optimal fits. Test optical predictions for future discrimination. Inspect existing visit records; fixed-opacity variable covering and fixed-temperature variable spots are affine in flux and cannot be distinguished solely by their functional temporal form. Use no false likelihood ranking or fitted confidence interval. Preserve native inputs, SHA256 hashes, source and tests.

Sampling amendment: fixed-cover varying-column model solves coverage analytically from g for each initial/final optical-depth pair, avoiding a coarse coverage grid that misses narrow solutions near the optical loss. Exact-g conditional slices, not full two-band tolerance volumes. Temporal extension uses the already saved 2018/2022 per-camera anchors and all other visit rows: fit column model to g/W1 anchors only, infer each subsequent optical depth from g alone, and evaluate W1 without refitting. Reject anchor residual >0.002 mag. Report model residual ranges; do not optimize models against evaluation visits.

Retrospective extension after temporalresults: test nominal5000K tau0=3 rescue against matched2018 absoluteg-W1 using nominalABg and savedVegaW1 zero points. No globalSEDfit or new atmosphere search;4500/5500K fixedgeometrysensitivity only. Calibrationassumption explicitly disclosed. This is not preregistered validation.
