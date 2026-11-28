# Q2 final relationship study outputs

Generated from the frozen Q2.1 six-curve output. This is descriptive relationship and effective-lag analysis; no prediction models are trained or evaluated.

- dose_response_curves.csv / dose_response_summary.csv: cumulative emission and gas-flow profiles over common cumulative-oxygen dose.
- gas_flow_by_oxygen_stage.csv: per-heat mean CO2/CO flow by four oxygen-dose stages.
- marginal_emission_yield_by_stage.csv: stage-wise incremental emission per oxygen volume.
- block_relationships.csv: non-overlapping 30-second increment associations.
- heat_emission_totals.csv: total emissions and emission/oxygen ratios by heat.
- lag_windows.csv / lag_summary.csv: window-level effective-lag estimates and summaries.
- lag_flow_bins.csv / lag_flow_slopes.csv: flow-conditioned lags and stage-adjusted within-heat slopes.
- within_heat_flow_contrasts.csv / _summary.csv: within-heat high-vs-low flow lag contrasts.
- synthetic_lag_recovery.csv: implementation sign/resolution check only.
- manifest.json: source fingerprint and analysis settings.

The analysis does not find robust evidence that higher lance flow yields shorter gas-signal delays. H0 is the main timebase; H1 is sensitivity only.
