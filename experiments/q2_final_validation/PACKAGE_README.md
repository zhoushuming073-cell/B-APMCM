# B-APMCM Q2.3 Final Validation and Freeze Package

This package completes the final, pre-specified validation of Q2.3. It does not introduce a fifth model-search round. Synthetic responses are used only to test the lag estimator; all contest-facing lag and flow values are computed from the ten official heats.

## One-command reproduction

Place the official `Data.zip` beside this README, or pass its path explicitly. Use a new, empty output directory:

```bash
python run_q23_final_validation.py --data-zip Data.zip --out q23_final_output --simulations 100 --seed 20261128
```

The command extracts the official workbook attachments, reconstructs the Q2.1 frozen oxygen/gas curves, runs the real-data Q2.3 analyses, executes 100 background resamples for each synthetic condition, audits archived results, and writes the final report. The output directory must be empty. Runtime depends on CPU; the 100-replicate experiment is split into ten deterministic batches of ten replicates to bound memory and execution time.

## Locked primary configuration

- Main time axis: H0; H1 is retained only as an existing alternative-time-axis sensitivity.
- Main window: 90 s, non-overlapping windows; sensitivity window: 180 s.
- Grid: 2 s; oxygen input uses the first two-thirds of each window; gas signal uses the full window.
- Main lag range: 0–30 s for 90-s windows; 0–60 s for 180-s windows.
- Trend removal: separate linear detrending of oxygen input and gas observation.
- Prewhitening: AR(1) coefficient estimated from detrended oxygen input and applied to both series.
- Peak: maximum absolute Pearson correlation; exact ties select the smallest nonnegative lag.
- Positive lag means oxygen input precedes gas observation.
- Random seed: 20261128; 100 replications per random condition.

## Outputs

`results/` contains:

- `actual_lag_windows.csv`, `actual_lag_summary.csv`, `actual_flow_bins.csv`: real ten-heat results.
- `matched_pairs_old_new.csv`, `matching_summary_old_new.csv`, `flow_progress_overlap_*.csv`: old/new process-progress alignment audit.
- `flow_lag_regression_old_new.csv`: original linear/quadratic stage-adjustment specifications under both progress definitions.
- `synthetic_recovery_windows.csv`, `synthetic_recovery_by_heat.csv`, `synthetic_recovery_equal_heat.csv`: known-delay recovery results.
- `noiseless_self_check_*.csv`: deterministic index/direction/boundary self-check.
- `zero_response_*.csv`, `zero_response_false_positive_summary.csv`: no-injected-response controls and calibration diagnostics.
- `round2_repaired_40_test_family.csv`: complete 40-test primary family and Holm correction; offset-shift sensitivities are retained separately.
- `*_archive_comparison.csv`, `historical_input_hash_audit.csv`: comparisons against the committed historical references.
- `figures/`: three core diagnostic plots.
- `Q2_3_FINAL_VALIDATION_AND_FREEZE.md`: final contest-facing conclusion and limitations.

The replicate-level synthetic CSV is included in the local output package. The GitHub repository should retain compact summaries and the deterministic code rather than duplicate the large raw replicate table.

## Interpretation boundary

The empirical residual blocks can contain genuine process variation; they are an empirical background, not pure measurement noise. The resulting semi-synthetic recovery experiment assesses algorithm behavior under constructed conditions. It is not an independent industrial transport-time measurement. The reported 14-s value is an absolute cross-correlation peak location (apparent time shift), not a validated physical transport delay.

## Final status

**封存状态：带明确限制封存。** The four previous research rounds are retained. The final validation establishes the limits of the estimator and the ten-heat flow-effect identification; it does not justify a monotone flow-delay model. No further Q2.3 model search is planned.
