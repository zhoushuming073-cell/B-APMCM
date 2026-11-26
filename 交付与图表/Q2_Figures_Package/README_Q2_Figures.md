# APMCM 2026 B题 — Question 2 Figure Package

This package contains figures generated from the ten-heat Question 2 experiment tables.

## Figure inventory

- `01_Heat_Process_Curves/`: one six-panel process-curve figure for each of the ten heats.
- `02_Cross_Heat_Comparison/`: ten-heat small-multiple comparisons of oxygen flow, cumulative oxygen, CO2/CO flow, and cumulative emissions.
- `03_Oxygen_Gas_Relationships/`: per-heat and pooled scatter plots of cumulative oxygen against gas flow and cumulative emissions, plus a per-heat Spearman heatmap.
- `04_Delay_Analysis/`: exploratory flow-bin lag, per-heat lag estimates, and lag-model MAE comparison.
- `05_Reconstruction_Validation/`: nested time-block validation comparisons for cumulative oxygen and gas-signal reconstruction.
- `06_Cumulative_Emissions/`: cumulative CO2/CO emissions by heat and oxygen/emissions comparisons.
- `07_Supporting_Data/`: CSV tables used for the figures and per-heat process-curve data.

## Important terminology

The available experiment tables contain `CO2_flow` and `CO_flow`, which are **instantaneous volumetric flow rates** (m3/s), not gas concentrations (e.g. vol%). Therefore, the figures label them as CO2/CO instantaneous volumetric flow. No concentration series is invented or inferred.

## Calculation and modeling notes

1. Time is taken from the supplied experiment curve tables (`time_s`), which use the official 0.5-second row interval. The special time-scale hypothesis for heat 6C02057 has not been resolved in this package; its plots retain the official 0.5-second time axis.
2. Cumulative oxygen is the recorded cumulative consumption (`oxygen_cumulative`, m3).
3. Instantaneous oxygen flow is an estimate derived from cumulative oxygen using the nonnegative integral-spline procedure (`oxygen_flow_spline`, m3/s). It is not an independently measured flow sensor series.
4. Cumulative CO2 and CO emissions are trapezoidal integrals of the recorded instantaneous volumetric flow series over the observed time window. No extrapolation beyond the recorded window is applied.
5. Scatter plots of cumulative oxygen against cumulative gas emissions can show high association partly because both quantities accumulate over time. They should not be interpreted as proof of a direct causal or linear physical law. The accompanying per-heat Spearman table includes both level and first-difference correlations.
6. Lag estimates are exploratory cross-correlation diagnostics. They are not independent physical transport-delay ground truth and do not establish causality.
7. Lag model MAE values compare offline alignment on a shared target/information setup; they are not equivalent to real-time forecasting accuracy.
8. Reconstruction validation metrics are outer time-block CV summaries. They assess interpolation/reconstruction of held-out observations, not agreement with independent physical sensors.

## Units

- Time: seconds (s)
- Oxygen consumption and cumulative emissions: cubic metres (m3)
- Instantaneous volumetric flow: cubic metres per second (m3/s)
- Lag: seconds (s)

## Suggested use in the contest paper

For the required process plots, use the ten individual six-panel figures or the cross-heat small multiples. For the relationship analysis, use the per-heat scatter plots and Spearman heatmap together. For the lag discussion, use the exploratory-lag figures together with the limitations stated above.
