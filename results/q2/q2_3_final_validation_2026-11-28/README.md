# Q2.3 Final Validation Results

These are compact repository-facing summaries. The full replicate-level simulation CSVs and one-command code are in the accompanying Q2_3_Final_Validation_and_Freeze.zip artifact.

- Actual ten-heat results: actual_lag_summary.csv, actual_flow_bins_h0_90.csv
- Flow identification: matching_summary_old_new.csv, flow_lag_regression_old_new.csv
- Estimator validation: synthetic_recovery_summary.csv, noiseless_self_check_summary.csv, zero_response_false_positive_summary.csv
- Reproducibility: historical_input_hash_audit.csv, curve_reconstruction_archive_comparison.csv, flow_profile_archive_comparison.csv, round2_archive_comparison.csv
- Historical 180-s table difference: prewhitened_version_comparison.csv

The interpretation is deliberately limited to apparent response time shifts. Do not call the 14-s peak a physical transport time.

Additional compact diagnostics: `synthetic_recovery_summary.csv` (100-replicate known-delay recovery), `signal_scale_diagnostics.csv` (actual/synthetic scale and ACF), `flow_progress_overlap_summary.csv`, `curve_equivalence_checks.csv`, and `simulation_batch_log.csv`.


## Q2.3 源码补归档（2026-11-28）

[完整复现源码与冻结输入](../../../experiments/q2_final_validation/) 已纳入仓库，包含原八个 Python 文件、依赖、配置及可校验恢复的 26 文件源包。不再依赖外部附件才能复跑。归档不重启模型研究，未重新执行 100 次模拟。 原始匹配表已找回：CO₂ 在 0.15 容差下的配对层中位差为 +10 s，炉次层中位差为 +8 s（炉次 Bootstrap 区间 [-26,20] s），旧精简 CSV 将配对层值误记为炉次层字段。以 [原始双统计量表](matching_summary_original_full.csv) 为准，历史文件保留；此前争议已由原归档直接澄清，无需重估，结论不变。
