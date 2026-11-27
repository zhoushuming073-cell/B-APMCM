# Q3 Physical Process Exploration (10 heats)

This directory records an exploratory, pre-modeling study of the ten official heats for Question 3.

## Research questions
1. Do the ten heats share recurring process stages, and can stage transitions be identified from the 2047-channel furnace-mouth spectrum?
2. Are CO and CO2 synchronous, lagged, or dynamically distinct?
3. Does the spectral-to-gas association vary with estimated oxygen-lance flow, cumulative oxygen progress, or process stage?

## Main findings
- CO2 crosses the defined onset threshold before CO in all 10 heats. Median lead: 77.5 s (range 24.0–172.5 s); exact two-sided sign test p=0.0020. The ordering persists at 5%, 20%, and 30% dynamic-range thresholds.
- A coarse three-phase description is directionally consistent for most heats, but BIC selects the maximum tested number of segments (K=6) in every heat. The three phases are descriptive, not a proven universal process ontology.
- Leave-one-heat-out balanced accuracy: 0.514 for 32-band spectrum only; 0.526 for full 2047-channel spectrum; 0.813 for time progress only; 0.743 for spectrum plus time. Current evidence does not support stable stage recognition from spectrum alone.
- CO and CO2 share long-run movement (median raw log-series zero-lag correlation 0.666), but their detrended zero-lag correlation is negative in all heats (median -0.194; stage 3 median -0.579). Cross-correlation peak lag is not stable across heats (-19.0 to +25.5 s).
- No stable flow-state interaction was detected. A small cumulative-oxygen-progress association for CO is borderline after Holm correction and is not robust to an alternative spectral PC1 measure.
- Heat 5 remains on the official 0.5 s time base in the primary analysis. A 1 s scenario is sensitivity analysis only.

## Methodological boundaries
- Instantaneous oxygen-lance flow is estimated by differencing cumulative oxygen; it is not a directly measured flow signal.
- Correlation and cross-correlation are descriptive, not causal evidence.
- Process-stage labels are constructed retrospectively using process variables; they must not be used as real-time predictors in a forecasting model.
- Heat 5/7 terminal gas collapses are retained in the primary analysis and only flagged in sensitivity checks.
- Full-spectrum stage classification is a diagnostic: per-heat unsupervised normalization and every-fifth-row sampling were used.

## Files
- `REPORT.md`: full methods, results, caveats, and implications.
- `summary_metrics.csv`: compact numerical summary.

The full exploratory archive (scripts, per-heat feature tables, detailed result tables, and figures) was generated separately. This repository record does not include the original official workbooks.
