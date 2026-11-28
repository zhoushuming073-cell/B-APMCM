# Integrated Q1-Q2-Q3 causal validation (2026-11-28)

## Scope
H0 is the official 0.5 s timebase. H1 is a sensitivity scenario (fifth heat at 1 s only). These are retrospective nested leave-one-heat-out experiments, not new independent furnace data. The 90 s lag analysis is Q2's pre-existing primary window, 180 s is its sensitivity analysis. The detected effective response lags are not verified physical transport delays. No Q1/Q2 frozen files have been changed.

## Results (macro R², 5 s feedback)

| Scenario / features | CO | CO2 |
|---|---:|---:|
| H0, Q2 90 s, integrated fusion | 0.964902 | 0.894178 |
| H0, Q2 90 s, old raw features, identical fold-specific lead | 0.963557 | 0.894062 |
| H0, Q2 180 s, integrated fusion | 0.907629 | 0.785609 |
| H1, Q2 90 s, integrated fusion | 0.975032 | 0.911497 |

H0 90 s leads 12–14 s, H0 180 s leads CO 32–38 s, CO2 36–40 s. The former improves forecast accuracy in part because the forecast is shorter; do not attribute the difference to the continuous reconstruction. The 90 s Q2 peak is weakly identified (0/71 BH-supported windows per gas), and the 180 s CO2 estimate also lacks strong support. Do not assert that the model forecasts true current bath chemistry with an exactly known delay. The model forecasts the future analyzer reading at an effective lead.

At H0 90 s, fusion macro R2 with 0/.5/5/15 s feedback: CO .846859/.969201/.964902/.953731; CO2 .696064/.899992/.894178/.878705. At 5 s, process+feedback .951391/.884196, initial-spectrum prior+feedback .963625/.916030, last-observation persistence .844419/.714566. Thus ongoing spectra provide modest benefit vs process+feedback and are not superior to the initial-spectrum prior for CO2. Fusion is a candidate for both gases, not a post-hoc gas-specific winner. Fifth heat H0 CO/CO2 .820133/.443686, H1 .922953/.619416 (different timebase/support). All ten heats included. No hardware validation of probe duty reduction.

## Method
Q1: log-intensity, one-sided shoulder-corrected temporal templates at 1/3/5 previous frames; fixed conservative agreement gate (not the original offline calibrated 99.9th-percentile detector); only 4 values repaired. All 14,805 saturated raw pixels retained as conditional lower bounds. Integrate the continuous piecewise-linear wavelength interpolant over 32 bands; time interpolation between already acquired frames follows Q1's bilinear rule. No unobserved future time interval is interpolated. This is an online adaptation, not an exact replay of Q1's bilateral repair or evidence of natural anomaly accuracy.

Q2: historical calibration uses the exact archived H0 per-heat spline parameters and two-iteration robust nonnegative integral cubic spline; Q2 k2 centered local means/PCHIP are used only for historical templates and historical lag estimation. For held-out runtime use the same cubic integral model class with 30 s trailing window, fixed 5 s spacing and penalty 10, NNLS and two IRLS steps. Endpoint cumulative estimates are made monotone via running maximum in the prediction features; these are not claimed to form a single globally mass-conserving spline. Raw gas remains the scoring target and sparse raw observations are the only online feedback. No future-centered test gas curve enters features.

Lag: reproduce Q2's 2 s grid, 90/180 s windows, linear detrending and AR(1) prewhitening with absolute Pearson peak over nonnegative lags. H0 window medians reproduce CO/CO2 14/14 and 32/36; H1 reproduces the same medians. Per-outer/inner fold, exclude the held-out heat and use median lag from training windows, rounded to nearest integer second. Window lag estimates computed independently per historical heat. The all-ten medians are descriptive only. For each fold/gas all compared families use the same lead and target times. No new monotonic flow-delay relationship imposed. H1 uses native 1 s grid only in E, and 5/15 s feedback intervals shared across timebases.

Q3: retained curve-transfer and ridge/heat-expert framework from previous expanded search, with a smaller fixed candidate set (alpha .03,.3,3; residual strengths .25,.5,1). Initial-frame-only optical prior; trailing runtime features; equal-heat ridge and expert fits; no PLS in this run. Feedback interval candidates 0/.5/5/15 s H0 and 0/5/15 s H1; gain grid 0–1.5 by .1. Model and gain chosen using inner LOHO on nine training heats. Compare full fusion, process-only, initial spectrum prior, ongoing spectrum-only, and gas persistence. No cross-heat external confirmation, no nested selection of evaluation horizon/feedback interval. Score detector timestamps >=60 s (coverage H0 12389/13589=91.17%; earlier detector times unscored). Feedback clock starts at detector 60 s and uses only sampled analyzer readings at or before publication. Early predictions run without feedback. Same 5 s candidate across comparisons.

## Independent audit
The verifier recomputes every metric from saved predictions and raw XLSB-derived gas arrays, checks exact target row sets and forecast lead, and validates every feedback clock. Raw input prefix tests at seven cuts per heat (70 H0 checks; additional 7 H1 fifth-heat oxygen checks) demonstrate prefix invariance. Model predictions re-executed on three truncated test feature arrays per fold (60 per run, 240 total). Event-driven sequential feedback is independently calculated for every selected model's scored origins and agrees exactly with vectorized output. Persistence timestamps are independently audited. Four runs: 1420 fold scores and 1,760,918 prediction records (repeated target records across methods/scenarios, not independent samples); zero timestamp violations; max metric difference <4e-15, model prefix difference <4e-14, feedback replay difference 0. The H0 90 s main run was repeated after timebase refactor and all 380 fold score rows were identical. Model prefix checks are sampled, not an independent per-time raw-model retraining. A prefix/clock test cannot establish unobserved external validity.

## Reproduction
Python 3.13, NumPy, pandas, SciPy, scikit-learn, pyxlsb, threadpoolctl. Linux process pool (fork). From this directory, using the verified official Data.zip (SHA256 cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727):

```bash
python extract.py --archive /path/to/Data.zip
OPENBLAS_NUM_THREADS=1 python prepare.py
OPENBLAS_NUM_THREADS=1 python prepare_h1.py
python prepare_raw_baseline.py
OPENBLAS_NUM_THREADS=1 python run.py --variant integrated --window 90
OPENBLAS_NUM_THREADS=1 python run.py --variant raw --window 90
OPENBLAS_NUM_THREADS=1 python run.py --variant integrated --window 180
OPENBLAS_NUM_THREADS=1 python run.py --scenario H1 --window 90
python verify.py
```

The zip includes the executable scripts, derived feature caches, lag-window tables, fold metrics, prediction records and QA. Raw XLSB-derived spectra are not duplicated; extract.py rebuilds them from the user's archive. Results can be regenerated without old experiment caches. Run directory names encode the scenario/window. The main run is `integrated_90`, alternative `integrated_180`, same-lead baseline `raw_90`, H1 `integrated_90_H1`. 
