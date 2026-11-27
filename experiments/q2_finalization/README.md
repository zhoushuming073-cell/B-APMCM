# Q2.1 finalization pipeline

This is the single reproduction entry point for the frozen Q2.1 outputs. It calls the versioned oxygen-flow V4 and gas-anomaly V1 scripts, reconstructs all six required curves for ten heats under H0 and H1, and writes compact validation tables.

From the repository root:

```bash
python experiments/q2_finalization/run.py
```

Default input: `data/interim/extracted/process.csv`.
Default output: `results/q2/final_freeze_2026-11-28/`.

The generated `final_curves.csv` contains 27,178 rows (13,589 source rows under each timebase scenario) and retains raw oxygen/CO2/CO observations alongside fitted curves. It is generated locally rather than versioned in Git because it is a multi-megabyte derived table. Compact summaries and QA tables are committed.

The pipeline asserts ten heats, 13,589 input rows, 27,178 final curve rows, 20 successful flow fits, nonnegative analytic flow, monotone fitted cumulative oxygen, mass-balance error below 0.1 m3, and no missing values in the final curve table.
