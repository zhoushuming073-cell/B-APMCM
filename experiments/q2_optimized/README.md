# Q2 legacy optimized experiment

This entry point is retained for historical reproduction of the earlier optimized run. It is **not** the canonical final Q2.1 six-curve pipeline.

For the frozen Q2.1 deliverable, run:

```bash
python experiments/q2_finalization/run.py
```

The final pipeline uses the nonnegative integral B-spline oxygen-flow estimate, the k2 CO2/CO reconstruction, cumulative integration, H0/H1 fifth-heat sensitivity, and the final anomaly review queue. See [`reports/final/Q2_1_FINAL_FREEZE.md`](../../reports/final/Q2_1_FINAL_FREEZE.md) and [`results/q2/final_freeze_2026-11-28/`](../../results/q2/final_freeze_2026-11-28/).

The historical command below remains available to reproduce the older experiment only:

```bash
python experiments/q2_optimized/run.py
```
