# Q3 causal gas-feedback audit

This branch records the V5 delayed-gas feedback reproduction and score-window audit.

## Rebuild the cache

Run from the repository root:

PYTHONPATH=src python experiments/q3_feedback_audit/prepare_cache.py --archive data/raw/official/Data.zip --out data/interim/q3_feedback_cache

The cache contains cumulative oxygen, 32 log-mean spectral bands, 32 line-contrast features, and CO/CO2 observations.

## Run the nested feedback sweep

python experiments/q3_feedback_audit/run_feedback_sweep.py --cache data/interim/q3_feedback_cache --out results/q3/feedback_audit/run_60s --score-start-s 60 --intervals 5,15,30,60 --gains 0.25,0.5,0.75,1.0

The outer protocol is leave-one-heat-out; model, lag, feedback interval, and gain selection occur within the remaining heats. The feedback schedule is causal: only an analyzer observation timestamped no later than the forecast publication time is used.

The audit is restricted to the frozen V5 q2_fixed lag-selection branch. It does not change the production Q3 entry point.
