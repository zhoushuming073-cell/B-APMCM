# Q2.3 Round 4 — Power and Identifiability Audit

This completed supplement evaluates whether the existing ten-heat design can detect a flow-dependent effective off-gas lag.

## Reproduce

From the repository root:

```bash
python experiments/q2_finalization/run.py
python experiments/q2_stage_matched_lag/v2/power_identifiability.py
```

Dependencies are the repository's existing NumPy, pandas and SciPy stack. No additional package is required.

## Scope

- Rebuilds the same 90 s / 180 s window lag estimates used in round three.
- Repeats within-heat 25/75 tail matching on cumulative oxygen progress.
- Compares linear and quadratic within-heat progress adjustment.
- Reports heat-clustered regression intervals and boundary-peak sensitivity.
- Runs 1,200 post-estimation simulations per effect cell, with two-window circular block residual resampling.

The simulation injects effects into estimated window-level lag outcomes. It is conditional power analysis, not a physical-process simulation or independent validation.

See `reports/final/Q2_3_STAGE_MATCHED_LAG_ROUND4_POWER_AUDIT.md` for interpretation and evidence boundaries.