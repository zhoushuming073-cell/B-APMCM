# Q2 process-gas anomaly review V1

## Scope

This detector reviews the ten heats' one-dimensional CO₂ and CO process-flow series. It extends the existing isolated one-point spike detector without deleting raw observations.

## Detection layers

1. **Isolated spike** — retains the existing nine-point flank comparison and raw-unit minimum threshold.
2. **Multi-point burst** — applies `log1p`, subtracts a 21-row rolling median, and searches for same-sign residual runs of 2–8 rows that return to baseline on both sides. The log-residual threshold is `max(3.5 × sigma, 0.018)`.
3. **Sustained level shift** — compares local linear extrapolations from 20 rows before and after a candidate boundary. It requires an absolute predicted level jump above `max(8 × sigma, 0.040)` and a sufficiently small pre/post slope mismatch. Nearby boundaries are clustered.

Here `sigma` is the 90th percentile of the absolute local log residual divided by 1.645, with a 0.003 floor. Thresholds are empirical screening thresholds, not calibrated probabilities.

## Current review queue

Across the 20 heat–gas series:

- 3 isolated-spike events (all CO)
- 10 multi-point bursts (8 CO, 2 CO₂)
- 38 sustained-level-shift boundaries (23 CO, 15 CO₂)
- 51 event candidates in total; 70 row-level flags

Every event is `pending` manual review. The new burst and level-shift detectors are flag-only. Raw observations are preserved and no new candidate is automatically corrected. The existing cleaned reconstruction path may still exclude legacy isolated-spike flags, but the original raw values remain available.

Two especially large late-series CO collapses are flagged:
- `6C02059`, row 1436: post-window median about 0.49 versus pre-window median about 34.24.
- `6C02057`, row 753: post-window median about 0.40 versus pre-window median about 31.83.

These may be genuine process stop/end-stage transitions rather than sensor faults. They must be checked against operating logs before any correction.

## Synthetic positive-control tests

Injected disturbances were added to the original series, using amplitudes of 1 or 4 gas-flow units:

| Disturbance | Duration | Amplitude | Event recall |
|---|---:|---:|---:|
| Isolated spike | 1 row | 1 | 78% |
| Isolated spike | 1 row | 4 | 99% |
| Multi-point burst | 3 rows | 1 | 73% |
| Multi-point burst | 3 rows | 4 | 84% |
| Multi-point burst | 5 rows | 1 | 72% |
| Multi-point burst | 5 rows | 4 | 69% |
| Sustained shift | 40+ rows | 1 | 50% |
| Sustained shift | 40+ rows | 4 | 100% |

Mean new flags outside the injected-event buffers were at most 0.4 per heat–gas test case. These are controlled positive-control results, not an independent false-positive estimate: the original series has no anomaly labels, and some baseline candidates may be real process transitions.

## Files

- `actual_review_queue.csv`: event-level manual review queue.
- `actual_row_flags.csv`: row-level flagged observations with raw values preserved.
- `detector_thresholds.csv`: per-series scale and threshold diagnostics.
- `injection_test_summary.csv`: event-recall summary by disturbance type.
- `injection_test_by_heat.csv`: per heat–gas injection results.
An example review figure was generated locally; raw values and event locations are fully represented in the CSV review queue.

## Timebase

Detection is performed in row space because all supplied records are evenly row-indexed. H0/H1 changes only the displayed time conversion; for the fifth heat, both time columns are retained in the review queue.

## Reproduction

From the repository root:

```bash
python experiments/q2_anomaly_review/v1/run.py \
  --input data/interim/extracted/process.csv \
  --out results/q2/anomaly_review/v1
```
