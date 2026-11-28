# Q2.2 lag-hypothesis robustness checks

This is the final supplementary verification for the Q2.2 hypothesis that higher oxygen-lance flow corresponds to a shorter effective off-gas response lag.

## Run

From the repository root, after regenerating the frozen Q2.1 curves:

```bash
python experiments/q2_finalization/run.py
python experiments/q2_lag_robustness/run.py
```

The script reads results/q2/final_freeze_2026-11-28/final_curves.csv and writes full actual and synthetic per-window tables, compact summaries, and a manifest.

## Checks

- 90 s and 180 s windows, with zero and half-window start offsets.
- Q2.1 analytic-derivative flow plus 30 s / 60 s centered differences of raw and fitted cumulative oxygen.
- Within-heat low/high flow contrasts and exact one-sided sign tests; Holm correction across 40 primary unshifted comparisons.
- Deterministic 20 s fixed-lag injection to diagnose estimator/window artifacts.

The synthetic test is an algorithm diagnostic, not a physical gas-process simulation. The study trains no prediction model and changes no source data.
