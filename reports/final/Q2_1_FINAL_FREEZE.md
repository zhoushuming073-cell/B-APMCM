# Q2.1 Final Freeze — 2026-11-28

> **单位语义补充审计（2026-11-28）：** 本报告冻结的平滑、拟合与数值重构仍可视为对所给 CO₂/CO 数值通道的曲线重构；但将其命名为 m³/s 流量并积分为绝对排放量，依赖题面单位解释正确。后续十炉条件氧平衡审计发现，若 CO/CO₂与氧枪累计供氧采用同一体积基准，COx 氧当量积分为氧枪供氧的 2.93–6.53 倍（中位数 5.77）。因此，绝对流量/排放解释暂列待核实，不覆盖原冻结数值结果。见 [Q3 测量语义与物理化学审计 V2](../working/Q3_PHYSICAL_CHEMISTRY_AUDIT_V2.md)。\n\n## Scope and source

Q2.1 delivers six time-dependent curves for each of ten heats: cumulative oxygen, effective instantaneous oxygen flow, CO2 flow, CO flow, cumulative CO2 emission, and cumulative CO emission.

The corrected source has 10 heats and 13,589 rows. Data.zip SHA-256: `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`.

- **H0 (main):** 0.5 s/row for all heats.
- **H1 (sensitivity):** 1.0 s/row for heat `6C02057` only; 0.5 s/row otherwise.

The fifth-heat timebase remains unresolved; H1 is not promoted to ground truth.

## Frozen models

### Oxygen flow

Fit a nonnegative cubic B-spline flow function and integrate it to the cumulative counter:

\[
q(t)=\sum_j c_jB_j(t),\quad c_j\geq0,\qquad
\widehat O(t)=O(0)+\int_0^tq(u)\,du.
\]

The instantaneous-flow curve is the analytic derivative, not a finite difference of the raw counter. Robust reweighting and a second-difference penalty are used. Knot spacing and penalty are selected by blocked within-heat validation. The one-standard-error envelope is parameter sensitivity, not a confidence interval.

Outer blocked validation of cumulative oxygen:

| Method | H0 MAE (m3) | H1 MAE (m3) |
|---|---:|---:|
| Nonnegative integral B-spline | **2.141** | **2.182** |
| Linear interpolation | 2.361 | 2.361 |
| PCHIP | 5.278 | 5.278 |

Across 20 heat/scenario fits, all fits succeeded; all flow values were nonnegative; fitted cumulative oxygen had no monotonicity violations. Maximum absolute trapezoidal mass-balance discrepancy was **0.0773 m3**, below the 0.1 m3 check threshold.

For heat `6C02057`, median effective flow is **8.979 m3/s (H0)** and **4.581 m3/s (H1)**. There is no independent high-frequency flowmeter truth, so the result is a physically constrained smoothed effective-flow estimate, not verified pointwise actual flow.

### CO2 and CO

The common primary model is **k2**: local means in a 2-second window at 2-second knot spacing, followed by shape-preserving PCHIP. Cumulative emissions use the exact antiderivative of the fitted PCHIP curve.

| Gas / scenario | k2 MAE (m3/s) | k2 RMSE (m3/s) |
|---|---:|---:|
| CO2, H0 | 0.04708 | 0.05970 |
| CO2, H1 | 0.04711 | 0.06003 |
| CO, H0 | 0.12829 | 0.16541 |
| CO, H1 | 0.13189 | 0.17145 |

CO2 k2 ranks first under both timebases. CO k4 is slightly better than k2 on the primary H1 outer-block MAE, but k2 wins the separate three-fold contiguous-gap robustness test under both H0 and H1. A single common k2 method is retained; k4 remains a sensitivity comparator.

Cumulative k2 emission totals differ from raw-flow trapezoidal totals by no more than **0.021%**. Three isolated CO candidates were found. Keeping or screening them leaves held-out MAE/RMSE unchanged to six decimals; full-curve sensitivity is at most 0.0848 m3/s and 0.2686 m3 cumulative. The final primary k2 curve retains all raw gas observations.

## Anomaly handling

| Candidate type | Events | Disposition |
|---|---:|---|
| Isolated one-point spike | 3 | Flag only |
| Multi-point burst | 10 | Flag only |
| Sustained level shift | 38 | Flag only |
| **Total** | **51** | **All pending manual review** |

These are candidates, not confirmed sensor faults. No source observation is automatically deleted or corrected. The detector's synthetic injection recall is 78%/99% for one-row spikes at amplitude 1/4; 73%/84% for three-row bursts; 72%/69% for five-row bursts; and 50%/100% for 40-row shifts. These are positive-control results, not real-world false-positive estimates.

Two major late CO collapses remain flagged: `6C02059` row 1436 (post-window median 0.49 vs pre-window 34.24) and `6C02057` row 753 (0.40 vs 31.83). They may be genuine process-end transitions and are not corrected without operating logs.

## Final QA

The final one-command pipeline generated 27,178 curve rows (13,589 under each timebase), 20 flow fits, 20/20 successful fits, 20/20 nonnegative checks, 20/20 cumulative-monotonicity checks, 20/20 mass-balance checks at 0.1 m3 tolerance, and 51 review events. Raw data changes: zero.

Run from repository root:

```bash
python experiments/q2_finalization/run.py
```

The full `final_curves.csv` is generated locally and not committed because it is a multi-megabyte derived table. Compact summaries, QA outputs, scripts, and this report are versioned.

## Freeze decision

**Q2.1 is frozen for the competition submission.** No further model search is planned. Remaining limitations are documented: fifth-heat timebase ambiguity, no independent instantaneous-flow ground truth, and no labelled anomaly truth.
