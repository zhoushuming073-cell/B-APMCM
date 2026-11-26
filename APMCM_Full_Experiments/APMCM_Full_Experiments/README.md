# APMCM Full-Data Experiments

## Input
- User-provided `/mnt/data/Data.zip`
- SHA-256: `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`
- Contains 10 XLSB heat workbooks.
- The raw files are not modified by these scripts.

## Main report
- `experiment_report.md`
- `data_audit.csv`
- `q1_full_reconstruction_summary.csv`
- `q1_gate_summary.csv`
- `q1_peak_masking_summary.csv`
- `q1_peak_masking_local_summary.csv`
- `q1_native_basis_summary.csv`
- `q1_anomaly_allheat.csv`
- `q1_anomaly_strata.csv`
- `q1_saturation_structure.csv`
- `q1_saturation_simulation_allheats.csv`
- `q2_flow_nested_cv.csv`
- `q2_flow_support_nested.csv`
- `q2_counter_updates_nested.csv`
- `q2_counter_simulation_allheats.csv`
- `q2_gas_nested_cv.csv`
- `q2_lag_by_heat.csv` (includes BH-adjusted q-values)
- `q2_lag_windows.csv`
- `q2_lag_loho_prediction.csv`
- `q2_phase_nested.csv`
- `q2_associations_nested.csv`
- `q2_cumulative_emissions.csv`

## Scripts
- `q1_anomaly_one.py`: per-heat synthetic anomaly screening.
- `q1_gate_one.py`: per-heat calibration/test of two-sided repair gate.
- `q1_shared_basis.py`: leave-one-heat-out peak masking with common-grid bases.
- `q1_native_basis.py`: leave-one-heat-out peak masking on native family grids.
- `q1_saturation_all.py`: actual saturation structure and simulated censoring tests.
- `q1_full_reconstruction.py`: provisional full-heat screening, gate, and saturation-template scenario.
- `q2_allheat_nested.py`: nested blocked CV for oxygen and gas reconstruction, plus process summaries.
- `q2_counter_simulation.py`: known-truth sample-and-hold simulation.
- `q2_delay_allheats.py`: per-heat lag diagnostics and nested LOHO prediction comparison.
- `make_figures_fast.py`: per-heat 3D spectral surfaces and six process curves.

## Generated figures and scenarios
- `figures/`: 10 spectral-surface PNGs and 10 six-panel process-curve PNGs.
- `cleaned_scenarios/`: per-heat NPZ files containing raw matrix, provisional scenario, score, screening mask, approved-repair mask, saturation mask, and saturation-repair mask.
- `*_candidate_top500.csv`: top 500 screened candidates by score per heat.

## Interpretation cautions
- Synthetic injections validate detection/recovery under the injected corruption model, not true anomaly precision.
- Peak masking uses algorithmically selected local maxima, not manually confirmed physical lines.
- The repair scenario is not a certified clean dataset; raw values and masks are retained.
- The lag analysis is exploratory and does not establish a causal or universal flow-delay function.
- The third question was not included in this experiment batch.

## Re-run order
1. Audit the supplied XLSB workbooks using `data_audit.csv`; `read_one.py` documents the per-workbook extraction routine.
2. Run Question 1 anomaly, masking, saturation, and scenario scripts.
3. Run `q2_allheat_nested.py`, `q2_counter_simulation.py`, `q2_delay_allheats.py`.
4. Run `make_figures_fast.py`.

Additional detailed exports:
- `q1_all_screen_candidates.csv`, `q1_approved_anomaly_repairs.csv`, `q1_saturation_repairs.csv`
- `q1_full_reconstruction_summary.csv`
- `q2_flow_lag_bins.csv`, `q2_loho_summary.csv`, `q2_time_third_macro.csv`, `q2_associations_macro.csv`
- `figures/all_heats_spectral_surface.pdf` and `figures/all_heats_process_curves.pdf`
- `亚太杯_完整数据实验更新.docx`: local V2.0 append-only research document. SharePoint replacement was blocked by a resource lock (HTTP 423), so the existing SharePoint file has not been overwritten.
