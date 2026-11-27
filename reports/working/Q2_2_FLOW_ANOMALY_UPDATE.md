> **Superseded for final submission by [`reports/final/Q2_1_FINAL_FREEZE.md`](../final/Q2_1_FINAL_FREEZE.md).** This working note records the V4 flow model and V1 anomaly detector. The final freeze report is authoritative for model choice, raw-value preservation, QA results, and remaining interpretation limits.

# Q2.2 update — instantaneous oxygen flow and process-gas anomaly review

Date: 2026-11-28

## 1. Instantaneous oxygen flow

### Model

The V4 reconstruction represents instantaneous oxygen flow as a nonnegative cubic B-spline and fits the cumulative oxygen counter through its analytic integral:

\[
q(t)=\sum_j c_jB_j(t),\quad c_j\geq0,\qquad
\widehat O(t)=O(0)+\int_0^tq(u)\,du.
\]

The flow estimate is therefore the analytic derivative of the fitted cumulative curve, not a finite difference of the raw counter. The coefficient fit uses robust reweighting and a second-difference smoothness penalty. A one-standard-error set of near-optimal hyperparameters provides a sensitivity envelope; this is not a confidence interval.

### Validation

Nested within-heat outer-block validation compares cumulative-oxygen reconstruction:

| Method | H0 mean MAE (m³) | H1 mean MAE (m³) |
|---|---:|---:|
| Linear interpolation | 2.361 | 2.361 |
| PCHIP | 5.278 | 5.278 |
| Nonnegative integral B-spline | **2.141** | **2.182** |

This validates short-gap reconstruction of cumulative oxygen, not the pointwise physical truth of instantaneous flow.

Across the 20 final heat/scenario fits:
- No negative flow estimates.
- No monotonicity violations in fitted cumulative oxygen.
- Maximum numerical mass-balance discrepancy between trapezoidal integration of q and the fitted cumulative increment: 0.0773 m³.

The median estimated flow is about 4.79–5.03 m³/s for nine heats. For heat `6C02057`, it is 8.98 m³/s under H0 and 4.58 m³/s under H1. The fifth-heat timebase remains unresolved; both scenarios are retained. The derivative is an effective smoothed estimate because no independent high-frequency flowmeter ground truth is available.

## 2. Process-gas anomaly review

The existing isolated-spike detector is retained. Two additional candidate detectors were added:
- Multi-point burst: same-sign log-residual run of 2–8 rows against a 21-row rolling-median baseline, with return-to-baseline checks.
- Sustained level shift: two-sided local linear extrapolation over 20 rows, requiring a sufficiently large level discontinuity and limited slope mismatch.

The current review queue contains 51 events: 3 isolated spikes, 10 multi-point bursts, and 38 sustained level-shift boundaries. They affect 70 rows. All remain pending manual review. New burst/shift candidates are flag-only; raw observations are preserved.

Two extreme late-series CO collapses deserve particular attention:
- `6C02059`, row 1436: post-window median 0.49 versus pre-window median 34.24.
- `6C02057`, row 753: post-window median 0.40 versus pre-window median 31.83.

These may be process stop/end-stage transitions, not sensor faults. Do not remove them without checking the operating log.

### Synthetic injection tests

| Disturbance | Duration | Amplitude | Event recall |
|---|---:|---:|---:|
| Isolated spike | 1 row | 1 unit | 78% |
| Isolated spike | 1 row | 4 units | 99% |
| Multi-point burst | 3 rows | 1 unit | 73% |
| Multi-point burst | 3 rows | 4 units | 84% |
| Multi-point burst | 5 rows | 1 unit | 72% |
| Multi-point burst | 5 rows | 4 units | 69% |
| Sustained shift | 40+ rows | 1 unit | 50% |
| Sustained shift | 40+ rows | 4 units | 100% |

These are positive-control tests on synthetic injections, not independent estimates of real-world false-positive rates. The source series has no labelled anomaly ground truth, and candidate level shifts can be genuine process transitions.

## 3. Decisions and remaining blockers

1. Use the analytic derivative as the primary effective-flow reconstruction; retain 10/30/60-second secant rates only as descriptive cross-checks.
2. Keep H0/H1 separate until the fifth heat's acquisition interval is verified.
3. Preserve the 51-event review queue. Only confirmed isolated spikes may follow the existing cleaned-reconstruction exclusion path; new burst and shift candidates must not be automatically deleted.
4. Check the two late-series collapses against process logs before deciding whether they are anomalies or shutdown events.
5. If the competition requires true instantaneous flow rather than a smoothed effective estimate, independent higher-frequency flowmeter data are needed.

## Reproduction

- `experiments/q2_instantaneous_flow/v4/run.py`
- `experiments/q2_anomaly_review/v1/run.py`

Both scripts read `data/interim/extracted/process.csv` by default and write versioned result folders under `results/q2/`.
