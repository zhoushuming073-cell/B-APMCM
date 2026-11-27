# Q2 prewhitened delay sensitivity experiment

This exploratory experiment checks whether local oxygen-flow-to-off-gas delay estimates survive detrending, AR(1) prewhitening, alternative dependence metrics, timebase sensitivity, and leave-one-heat-out evaluation.

## Run

Run from the repository root:

```powershell
python experiments/q2_prewhitened_delay/run.py
```

The script reads the ten tracked files under `data/processed/process_curves/`. By default, it writes four tables to `results/q2/prewhitened_delay/rerun/`:

- `window_lags.csv`: one row per heat, gas, window, timebase, and dependence metric
- `summary.csv`: multiplicity-adjusted support, flow-lag slopes, heat-cluster bootstrap intervals, and LOHO scores
- `loho_scores.csv`: per-held-out-heat fixed-mean versus flow-linear MAE
- `flow_bins.csv`: descriptive low/middle/high flow-bin lag summaries

To choose another output directory:

```powershell
python experiments/q2_prewhitened_delay/run.py --output results/q2/prewhitened_delay/rerun
```

The committed files under `results/q2/prewhitened_delay/tables/` are archived experiment artifacts. A fresh run writes to `rerun/` and does not overwrite them. The Python port has not yet been independently validated end to end, so a successful execution should not be described as an independent reproduction until its outputs have been compared with the archived tables.

This experiment is intentionally not included in the default `scripts/run_all.py` batch. Add it to the batch only after independent reproduction and review.

## Protocol

- H0 retains the official 0.5 s per row clock.
- H1 is a sensitivity case only: heat 6C02057 uses 1.0 s per row, so its derived oxygen flow is halved.
- Signals are interpolated to a 2 s grid.
- Non-overlapping 90 s and 180 s windows are used. The input-flow segment occupies the first two-thirds of each window; the gas segment spans the full window.
- Both signals are locally linearly detrended. An AR(1) coefficient is estimated from the input-flow residual and the same filter is applied to the gas residual.
- Pearson, Spearman, and distance correlation are each scanned over 0–30 s (90 s windows) or 0–60 s (180 s windows).
- A circular-shift max-statistic supplies a local null distribution. Benjamini–Hochberg correction is applied separately within each metric × timebase × window × gas family.
- Flow-lag slopes are estimated after centering within heat, with equal heat weighting. Confidence intervals resample heats, not individual time points. The predictive check is leave-one-heat-out.

## Interpretation boundary

This is a sensitivity experiment, not a replacement for the frozen Q2 results. Circular-shift p-values have coarse resolution in short windows, and ten heats—not the thousands of within-heat observations—are the independent experimental units. A non-significant result does not prove that no physical response exists; it means this dataset does not identify a stable flow-conditioned lag with the tested protocol.
