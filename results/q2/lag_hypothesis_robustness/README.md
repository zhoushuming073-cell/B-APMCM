# Q2.2 lag-hypothesis robustness outputs

Committed compact evidence:
- key_actual_results.csv: selected actual CO/CO2 lag results for the frozen analytic flow, shifted windows, and 180 s multi-estimator sensitivity.
- key_synthetic_results.csv: deterministic 20 s constant-lag injection results.
- manifest.json: input fingerprint, heat/scenario coverage and test settings.

The experiment script regenerates full per-window tables (actual_lag_windows.csv and synthetic_fixed_lag_windows.csv) and complete summary tables locally. They are not committed to keep the repository compact.

The main finding is that the apparent inverse flow-lag signal is sensitive to window origin and flow estimator. No primary comparison survives Holm correction. Q2.2 is formally archived; no further model search is planned.
