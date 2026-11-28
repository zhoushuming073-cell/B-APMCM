# Q2.3 Final Validation Results

These are compact repository-facing summaries. The full replicate-level simulation CSVs and one-command code are in the accompanying Q2_3_Final_Validation_and_Freeze.zip artifact.

- Actual ten-heat results: actual_lag_summary.csv, actual_flow_bins_h0_90.csv
- Flow identification: matching_summary_old_new.csv, flow_lag_regression_old_new.csv
- Estimator validation: synthetic_recovery_summary.csv, noiseless_self_check_summary.csv, zero_response_false_positive_summary.csv
- Reproducibility: historical_input_hash_audit.csv, curve_reconstruction_archive_comparison.csv, flow_profile_archive_comparison.csv, round2_archive_comparison.csv
- Historical 180-s table difference: prewhitened_version_comparison.csv

The interpretation is deliberately limited to apparent response time shifts. Do not call the 14-s peak a physical transport time.

Additional compact diagnostics: `synthetic_recovery_summary.csv` (100-replicate known-delay recovery), `signal_scale_diagnostics.csv` (actual/synthetic scale and ACF), `flow_progress_overlap_summary.csv`, `curve_equivalence_checks.csv`, and `simulation_batch_log.csv`.
