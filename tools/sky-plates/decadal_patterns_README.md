# J0754 decadal pattern checks

`decadal_patterns.py` reads existing saved J0754 photometry and produces annual summaries, same-series local event contrasts, matched-control comparisons, actual-cadence historical waveform replay and a bounded retrospective 14.8-year schedule check. It does not download or modify source data. Run with the existing runtime/sky-plates/venv Python from this directory. Outputs are under results/sky-plates-j0754-decadal-patterns-20261109; the report, protocol and JSON document exact inputs and limitations.

Annual minima use actual median observation times. The recurrence schedule is constructed from selected endpoints, not estimated with a calibrated statistical search. Historical gaps are not bright states. Optional linear-drift replay fits are essential sensitivity results and must accompany offset-only comparisons. Different camera/filter baselines are not physically merged.

Tests: python -m unittest test_decadal_patterns. Full component regression suite: python -m unittest discover -s . -p "test_*.py".
