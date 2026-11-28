# Q2.3 Round 4 results

This directory stores compact, inspectable outputs from the fourth post-freeze audit.

| File | Contents |
|---|---|
| `matching_summary.csv` | Third-round 25/75 tail matching reproduced under H0/H1 |
| `matched_pairs.csv` | Individual high/low window pairs and progress gaps |
| `heat_slope_summary.csv` | Heat-level continuous partial slopes and bootstrap summaries |
| `pooled_cluster_regression.csv` | Linear/quadratic stage-adjusted coefficients with heat-clustered intervals; all/interior windows |
| `power_curve_linear_stage.csv` | Conditional power under linear stage adjustment |
| `power_curve_quadratic_stage.csv` | Conditional power under quadratic stage adjustment |
| `mde_summary.csv` | First tested effect size reaching 80% power |
| `design_diagnostics.csv` | Per-heat flow/progress relationship and residual-flow variation |
| `heat_diagnostics.csv` | Per-heat partial slope diagnostics |
| `input_reconstruction_qa.json` | Source fingerprint, curve hash comparison and numerical reproduction checks |
| `manifest.json` | Experiment protocol, simulation settings and evidence limitations |

The complete 856-row window table is generated locally as `lag_windows.csv` by the experiment script and is not versioned. Reproduce from the repository root after generating Q2.1 curves:

```bash
python experiments/q2_finalization/run.py
python experiments/q2_stage_matched_lag/v2/power_identifiability.py
```
