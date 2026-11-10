# Broader J0754 dust sensitivities

Run with the existing sky-plates Python environment:

python broader_dust.py --root <Lumen root> --out <results directory>

Requires existing CK04/SVO model inputs, Draine tables, original model comparison and chromaticity results. No downloads. Results compare nonzero-baseline covering/column models, gray/selective and smooth power-law extinction, and dust-plus-spots. g/W1 select parameters; W2 is retrospective check only. Grid envelopes anchor g exactly and allow W1 +/-0.02mag, not full confidence regions. The temporal test fits only two existing anchors and evaluates other W1 visits from g. Reuses an existing scientific component, not a new project. See output protocol/report for parameter ranges and limits.
