# Q2.2 final relationship study

This workflow studies the relationships requested in Problem B Question 2. It is descriptive and lag-characterisation research only: **it does not train, tune, or evaluate prediction models**.

## Run

From the repository root:

    python experiments/q2_finalization/run.py
    python experiments/q2_relationship_study/run.py

The second command reads the frozen Q2.1 final_curves.csv and writes compact and full-resolution CSV tables to results/q2/relationship_study_final/.

## Scope

- CO2/CO instantaneous gas-flow profiles over common cumulative-oxygen dose.
- Cumulative CO2/CO emissions as functions of cumulative-oxygen dose.
- Stage-wise marginal emission per oxygen volume.
- Non-overlapping block associations using increments, not raw monotone levels.
- Effective lag estimates under low/middle/high oxygen-flow conditions, with H0/H1 timebase sensitivity and within-heat stage adjustment.

No source values are modified. H0 is the main timebase; H1 is only a fifth-heat sensitivity case. The lag is an effective process-and-measurement response delay, not pure pipe transit time.

## Important interpretation

The analysis does not find robust evidence that higher lance flow produces shorter gas-signal delays. Report the measured flow-stratified lags and this limitation honestly; do not force a monotone delay law to match the problem statement.
