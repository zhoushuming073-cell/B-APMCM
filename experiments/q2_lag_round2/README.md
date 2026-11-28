# Q2.3 Round 2 — Unified lag search and estimator sensitivity

This isolated branch adds a second-round sensitivity experiment for the effective response lag between oxygen-lance flow and CO/CO2 off-gas flow. It does not modify the frozen Q2.1 or Q2.2 results and does not involve Question 3.

## Design

- Common nonnegative lag search: 0–30 s at 2 s increments for 90, 120 and 180 s windows.
- Each window uses the first approximately two-thirds of the oxygen-flow series as input and the full fixed gas-response window as output.
- Compare positive-only and absolute correlation peaks.
- Compare Pearson and Spearman association.
- Compare detrending alone with AR(1) input prewhitening (the same fitted filter is applied to the response).
- Compare Q2.1's analytic derivative with 2, 30 and 60 s centered finite differences of raw and fitted cumulative oxygen.
- H0 is primary; H1 is the fifth-heat timebase sensitivity scenario.
- Flow-stratified contrasts are computed within heat. Bootstrap resampling uses heat, not windows, as the independent unit.

## Run

From the repository root:

    python experiments/q2_finalization/run.py --input data/interim/extracted/process.csv --out results/q2/final_freeze_2026-11-28
    python experiments/q2_lag_round2/run.py --input results/q2/final_freeze_2026-11-28/final_curves.csv --output results/q2/lag_round2_2026-11-28

## Interpretation

The pre-specified comparable primary configuration is H0, 90 s, analytic derivative, Pearson, AR(1) prewhitening, absolute peak. All other combinations are sensitivity analyses. A peak at 30 s is right-censored by the common search range and must not be interpreted as an exact 30 s delay.

The experiment is descriptive. It does not identify pure pipe transit time and does not presuppose that higher oxygen flow must produce shorter lag.
