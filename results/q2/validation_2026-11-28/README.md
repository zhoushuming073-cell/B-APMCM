# Q2.1 Reconstruction Validation — 2026-11-28

This directory contains the machine-readable result tables from the Q2.1 validation rerun.

## Data and timebase scenarios

- Source: ten corrected XLSB attachments supplied in `Data.zip`; 13,589 time points in total.
- Source archive SHA-256: `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`.
- **H0:** 0.5 seconds per row for all ten heats.
- **H1:** 1.0 second per row for the fifth heat (`6C02057`) only; 0.5 seconds per row for the other nine heats.

## Gas reconstruction models

- `linear`: raw observations with linear interpolation.
- `clean_linear`: isolated-spike flagging followed by linear interpolation.
- `pchip`: isolated-spike flagging followed by shape-preserving PCHIP interpolation.
- `k2`, `k4`, `k8`: isolated-spike flagging, time-bin averaging over 2, 4, or 8 seconds, followed by PCHIP interpolation.
- `smoothing_spline`: cubic smoothing spline; smoothing parameter selected using inner validation blocks only.

`gas_model_ranking.csv` reports the outer leave-one-heat-out validation results. Each heat contributes 32 held-out points from two contiguous blocks. `gas_gap_cv_ranking.csv` reports a separate three-fold contiguous-gap stress validation. MAE and RMSE are in m³/s. Lower is better.

## Other validation tables

- `gas_peak_reconstruction_check.csv`: held-out local peak magnitude and peak-location errors.
- `emission_totals.csv`: cumulative gas integrals from each full-series reconstruction, compared with the raw-data trapezoidal integral.
- `oxygen_flow_ranking.csv`, `oxygen_flow_per_heat.csv`, `oxygen_flow_full_refit.csv`: cumulative-oxygen holdout validation and nonnegative integral B-spline flow refits. The cumulative-oxygen error does **not** independently validate true instantaneous oxygen flow.
- `anomaly_*`: observed candidate points and synthetic anomaly-injection tests. The current detector is designed for isolated one-point spikes and should not be treated as a general regime-shift detector.
- `synthetic_flow_stress_summary.csv`: synthetic step, pulse, and oscillatory flow tests using empirical counter-update patterns; these are stress tests, not observed ground truth.
- `data_quality_audit.csv`: basic missingness, range, and monotonicity checks.

## Interpretation limits

There is no independent ground truth for gas values between the original observations or for instantaneous oxygen flow. Model rankings therefore describe interpolation/reconstruction performance under the specified held-out designs, not sensor-level physical truth. H0/H1 are retained as separate sensitivity scenarios; the fifth-heat timebase is not declared resolved.

The full 20-sheet workbook is distributed separately as `Q21_Validation_Results.xlsx`; this GitHub directory stores the CSV tables and methodological note for version control and diffability.