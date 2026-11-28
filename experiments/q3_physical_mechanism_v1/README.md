# Q3 Physical Mechanism V1

This experiment investigates three hypotheses using the ten official heats:

1. Whether the terminal CO2 rise / CO decline is consistent with a falling total carbon-bearing gas flow and a rising apparent CO2 share.
2. Whether local spectral-to-gas alignment differs by retrospective process stage.
3. Whether physically localized spectral features and causal momentum features improve leave-one-heat-out prediction, especially around the terminal reversal.

## Data and timing

- Input: official corrected Data.zip, SHA-256 cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727.
- Main timebase: 0.5 s per row for all heats (H0).
- Ten heats, 13,589 rows, 2,047 native wavelength channels per row.
- CO2 and CO columns are treated as instantaneous volume flow rates (m3/s), not ppm concentrations.
- Raw files are not modified. Five-second median filtering is used only for diagnostic curves and is applied along the time axis only (11 x 1), never across the CO2 and CO columns.

## Key methods

- Tail window: final 60 s before the first sustained terminal collapse in heats 6C02057 and 6C02059; otherwise the last 60 s of the record. Sensitivity windows: 30, 45, 60 and 90 s.
- Apparent CO2 share: CO2 / (CO2 + CO). This is not a full industrial post-combustion ratio.
- Four-segment gas-state fit: a diagnostic segmentation of log(total gas flow) and logit(apparent CO2 share), with four segments fixed to test the team's stage hypothesis. Labels are retrospective and are not used as online predictors.
- Stage-lag diagnostic: 5 s differences of spectral brightness and gas flow, linearly detrended within retrospective stages, searched over 0–60 s in 5 s increments. These are statistical offsets, not pipe transport times.
- K I feature: local doublet contrast around 766.490 and 769.896 nm relative to nearby continuum windows. It is a candidate state feature, not a potassium concentration measurement.
- Prediction ablation: nested leave-one-heat-out Ridge comparisons; inner folds choose regularization, and training heats choose the spectral lead. Results are exploratory and not directly comparable with archived V5 scores.

## Results

See the research report at reports/working/Q3_PHYSICAL_MECHANISM_V1.md.

Frozen compact outputs:
- tail_dynamics.csv
- tail_window_sensitivity.csv
- stage_summary.csv
- stage_lag_summary.csv
- spectral_event_features.csv
- ablation_macro.csv
- ablation_fold_delta.csv
- decomposition_macro.csv
- tail_prediction_slopes.csv

## Limitations

- Ten heats have been repeatedly explored; leave-one-heat-out is not independent confirmation.
- No direct bath carbon/silicon/manganese, slag FeO, mouth O2, lance height, temperature, pipe flow or analyzer step-response data are available.
- Stage labels are gas-defined retrospective diagnostics, not deployable state labels.
- The first notebook implementation has not yet been converted into a clean command-line reproduction script. Do not treat this README as a claim that the full run is reproducible from a single command. The next revision must archive and execute that script before model promotion.
