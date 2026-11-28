#!/usr/bin/env python3
"""Q2.3 Round 4: design-conditional power and identifiability audit.

Run from the repository root after regenerating the frozen Q2.1 curves:
    python experiments/q2_finalization/run.py
    python experiments/q2_stage_matched_lag/v2/power_identifiability.py

The simulation is post-estimation: it injects known effects into window-level lag
estimates while preserving empirical residual shapes. It is not a physical-process
simulation and does not establish a causal effect.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest, t, ttest_1samp

HEATS = [
    "6B03205", "6B03206", "6B03270", "6B03271", "6C02057",
    "6C02058", "6C02059", "6C02060", "6C02107", "6C02108",
]
SCENARIOS = ["H0", "H1"]
GASES = ["CO2", "CO"]
WINDOWS = [90, 180]
METRICS = ["pearson", "spearman"]
EFFECTS = [0, -2, -4, -6, -8, -10, -12, -16]
N_SIM = 1200
SEED = 20261128


def corr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    x = x - x.mean()
    y = y - y.mean()
    den = np.sqrt((x @ x) * (y @ y))
    return float(x @ y / den) if den > 0 else np.nan


def assoc(x, y, metric):
    if metric == "pearson":
        return corr(x, y)
    return corr(pd.Series(x).rank(method="average"), pd.Series(y).rank(method="average"))


def detrend(x):
    x = np.asarray(x, float)
    tt = np.arange(len(x), dtype=float)
    X = np.column_stack([np.ones(len(x)), tt])
    return x - X @ np.linalg.lstsq(X, x, rcond=None)[0]


def whiten(x, y):
    x = detrend(x)
    y = detrend(y)
    den = x[:-1] @ x[:-1]
    phi = np.clip((x[:-1] @ x[1:]) / den if den > 0 else 0.0, -0.98, 0.98)
    return x[1:] - phi * x[:-1], y[1:] - phi * y[:-1], float(phi)


def estimate_windows(curves):
    rows = []
    for metric in METRICS:
        for scenario in SCENARIOS:
            for window in WINDOWS:
                n_input = window // 3
                maxlag = window // 6
                for heat in HEATS:
                    f = curves[(curves.scenario == scenario) & (curves.heat == heat)].sort_values("time_s")
                    if f.empty:
                        raise ValueError(f"Missing curve: {scenario}/{heat}")
                    tt = f.time_s.to_numpy(float)
                    q = f.oxygen_flow_analytic_m3_s.to_numpy(float)
                    O = f.oxygen_cumulative_fit_m3.to_numpy(float)
                    for start in np.arange(0, tt[-1] - window + 1e-8, window):
                        grid = start + 2 * np.arange(window // 2)
                        qx = np.interp(grid[:n_input], tt, q)
                        progress = (np.interp(start + window / 2, tt, O) - O[0]) / (O[-1] - O[0])
                        for gas in GASES:
                            gy = np.interp(grid, tt, f[f"{gas}_flow_k2_m3_s"].to_numpy(float))
                            if np.std(qx) < 1e-7 or np.std(gy) < 1e-7:
                                continue
                            x, y, phi = whiten(qx, gy)
                            if np.std(x) < 1e-9 or np.std(y) < 1e-9:
                                continue
                            values = [abs(assoc(x, y[k:k + len(x)], metric)) for k in range(maxlag + 1)]
                            k = int(np.nanargmax(values))
                            rows.append({
                                "metric": metric, "scenario": scenario, "window_s": window,
                                "heat": heat, "gas": gas, "start_s": float(start),
                                "flow_m3_s": float(qx.mean()), "oxygen_progress": float(progress),
                                "lag_s": float(2 * k), "peak_abs_assoc": float(values[k]),
                                "signed_assoc_at_peak": float(assoc(x, y[k:k + len(x)], metric)),
                                "ar1_phi_from_oxygen": phi, "boundary_peak": bool(k == maxlag),
                            })
    return pd.DataFrame(rows)


def stage_matrix(progress, degree):
    p = np.asarray(progress, float)
    return np.column_stack([p ** k for k in range(degree + 1)])


def residualize(y, X):
    beta = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)[0]
    return np.asarray(y, float) - X @ beta, beta


def exact_sign_p(values):
    values = np.asarray(values, float)
    nz = values[np.abs(values) > 1e-10]
    if len(nz) == 0:
        return 1.0
    return float(binomtest(int((nz < 0).sum()), len(nz), 0.5, alternative="greater").pvalue)


def round3_pairs(g, caliper=0.15):
    """Reproduce round-three strict-tail, one-to-one progress matching."""
    out = []
    for heat, h0 in g.groupby("heat", sort=False):
        h = h0.sort_values("start_s").reset_index(drop=True)
        qlo, qhi = np.quantile(h.flow_m3_s, [0.25, 0.75])
        lo = h[h.flow_m3_s < qlo].sort_values("oxygen_progress").reset_index()
        hi = h[h.flow_m3_s > qhi].sort_values("oxygen_progress").reset_index()
        if lo.empty or hi.empty:
            continue
        dist = np.abs(hi.oxygen_progress.to_numpy()[:, None] - lo.oxygen_progress.to_numpy()[None, :])
        nh, nl = dist.shape
        cost = np.full((nh, nl + nh), caliper + 1e-8)
        cost[:, :nl] = np.where(dist <= caliper, dist, caliper + 1.0)
        for i in range(nh):
            cost[i, nl + i] = caliper + 1e-8
        rr, cc = linear_sum_assignment(cost)
        pairs = []
        for i, j in zip(rr, cc):
            if j >= nl or dist[i, j] > caliper:
                continue
            pairs.append({
                "heat": heat, "high_index": int(hi.iloc[i]["index"]),
                "low_index": int(lo.iloc[j]["index"]),
                "progress_gap": float(dist[i, j]),
                "flow_diff_m3_s": float(hi.iloc[i].flow_m3_s - lo.iloc[j].flow_m3_s),
                "lag_diff_high_minus_low_s": float(hi.iloc[i].lag_s - lo.iloc[j].lag_s),
            })
        if pairs:
            out.append((heat, h, pairs))
    return out


def pooled_cluster_regression(g, degree, boundary_mode):
    """Heat-specific stage polynomial + standardized within-heat flow residual.

    Returns CR1 cluster-robust coefficient and t interval with heat-level df.
    This implementation uses NumPy only, avoiding an additional statsmodels dependency.
    """
    d = g.copy()
    if boundary_mode == "interior":
        d = d[~d.boundary_peak].copy()
    d = d.sort_values(["heat", "start_s"]).reset_index(drop=True)
    n = len(d)
    if n < 20:
        return None
    heat_order = list(d.heat.drop_duplicates())
    stage_cols = []
    xresid = np.zeros(n)
    for heat in heat_order:
        ix = np.flatnonzero(d.heat.to_numpy() == heat)
        h = d.iloc[ix]
        Xp = stage_matrix(h.oxygen_progress.to_numpy(), degree)
        xr, _ = residualize(h.flow_m3_s.to_numpy(), Xp)
        sd = float(np.std(xr, ddof=1))
        if sd > 1e-10:
            xresid[ix] = xr / sd
        for k in range(degree + 1):
            col = np.zeros(n)
            col[ix] = h.oxygen_progress.to_numpy(float) ** k
            stage_cols.append(col)
    X = np.column_stack(stage_cols + [xresid])
    y = d.lag_s.to_numpy(float)
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    resid = y - X @ beta
    meat = np.zeros((X.shape[1], X.shape[1]))
    groups = d.heat.to_numpy()
    for heat in heat_order:
        ix = np.flatnonzero(groups == heat)
        score = X[ix].T @ resid[ix]
        meat += np.outer(score, score)
    G = len(heat_order)
    K = X.shape[1]
    correction = (G / (G - 1)) * ((n - 1) / max(n - K, 1))
    cov = correction * xtx_inv @ meat @ xtx_inv
    se = float(np.sqrt(max(cov[-1, -1], 0.0)))
    b = float(beta[-1])
    df = G - 1
    crit = float(t.ppf(0.975, df))
    p_less = float(t.cdf(b / se, df)) if se > 0 else 1.0
    return {
        "degree": degree, "boundary_mode": boundary_mode, "n_windows": n,
        "n_heats": G, "beta_s_per_within_heat_sd": b, "cluster_se": se,
        "ci_low": b - crit * se, "ci_high": b + crit * se,
        "one_sided_p_beta_lt_0": p_less,
    }


def heat_slopes(g, degree):
    rows = []
    for heat, h0 in g.groupby("heat", sort=False):
        h = h0.sort_values("start_s").reset_index(drop=True)
        Xp = stage_matrix(h.oxygen_progress.to_numpy(), degree)
        if len(h) <= degree + 2:
            continue
        x, _ = residualize(h.flow_m3_s.to_numpy(), Xp)
        sx = float(np.std(x, ddof=1))
        if sx < 1e-10:
            continue
        x = x / sx
        y, _ = residualize(h.lag_s.to_numpy(), Xp)
        slope = float((x @ y) / (x @ x))
        rows.append({"heat": heat, "n_windows": len(h), "slope_s_per_sd": slope,
                     "flow_resid_sd_m3_s": sx, "boundary_rate": float(h.boundary_peak.mean())})
    return pd.DataFrame(rows)


def circular_block_boot(e, rng, block=2):
    n = len(e)
    out = []
    while len(out) < n:
        s = int(rng.integers(0, n))
        out.extend([e[(s + k) % n] for k in range(block)])
    return np.asarray(out[:n], float)


def prepare_power_groups(g, degree):
    groups = []
    for heat, h0 in g.groupby("heat", sort=False):
        h = h0.sort_values("start_s").reset_index(drop=True)
        if len(h) <= degree + 2:
            continue
        Xp = stage_matrix(h.oxygen_progress.to_numpy(), degree)
        xraw, _ = residualize(h.flow_m3_s.to_numpy(), Xp)
        sx = float(np.std(xraw, ddof=1))
        if sx < 1e-10:
            continue
        x = xraw / sx
        stage, _ = residualize(h.lag_s.to_numpy(), Xp)
        stage = h.lag_s.to_numpy() - stage  # fitted stage component
        Xfull = np.column_stack([Xp, x])
        e, _ = residualize(h.lag_s.to_numpy(), Xfull)
        groups.append({"heat": heat, "h": h, "x": x, "stage": stage, "e": e})
    return groups


def power_simulation(g, degree, matched, n_sim, seed):
    rng = np.random.default_rng(seed)
    groups = prepare_power_groups(g, degree)
    pair_map = {heat: pairs for heat, _, pairs in matched}
    for z in groups:
        z["pairs"] = pair_map.get(z["heat"], [])
        # Pre-generate block-bootstrap residual paths once; common random numbers
        # make power differences across effect sizes less noisy.
        h, e = z["h"], z["e"]
        n = len(h)
        n_blocks = int(np.ceil(n / 2))
        starts = rng.integers(0, n, size=(n_sim, n_blocks))
        offsets = np.arange(2)[None, None, :]
        ix = (starts[:, :, None] + offsets) % n
        ix = ix.reshape(n_sim, -1)[:, :n]
        z["boot_e"] = e[ix]
    rows = []
    n_cont = len(groups)
    n_match = sum(bool(z["pairs"]) for z in groups)
    def sign_reject(matrix):
        # Exact one-sided binomial sign test. The rejection threshold is
        # precomputed for this fixed number of heat clusters.
        if matrix.shape[1] == 0:
            return np.zeros(matrix.shape[0], dtype=bool)
        n = matrix.shape[1]
        k = np.sum(matrix < 0, axis=1)
        critical = next((j for j in range(n + 1)
                         if binomtest(j, n, .5, alternative="greater").pvalue <= .05), n + 1)
        return k >= critical
    for effect in EFFECTS:
        slope_cols = []
        match_cols = []
        for z in groups:
            h, x, stage, boot_e = z["h"], z["x"], z["stage"], z["boot_e"]
            ysim = np.clip(np.round((stage[None, :] + effect * x[None, :] + boot_e) / 2) * 2, 0, 30)
            bh = (ysim - stage[None, :]) @ x / (x @ x)
            slope_cols.append(bh)
            if z["pairs"]:
                dd = np.column_stack([ysim[:, p["high_index"]] - ysim[:, p["low_index"]]
                                      for p in z["pairs"]])
                match_cols.append(np.median(dd, axis=1))
        slope_mat = np.column_stack(slope_cols) if slope_cols else np.empty((n_sim, 0))
        match_mat = np.column_stack(match_cols) if match_cols else np.empty((n_sim, 0))
        reject_sign = sign_reject(slope_mat)
        if n_cont > 1:
            means = slope_mat.mean(axis=1)
            sds = slope_mat.std(axis=1, ddof=1)
            tstat = means / (sds / np.sqrt(n_cont))
            crit = t.ppf(.05, n_cont - 1)
            reject_t = tstat < crit
        else:
            reject_t = np.zeros(n_sim, dtype=bool)
        reject_match = sign_reject(match_mat)
        rows.append({
            "effect_s_per_within_heat_sd": effect, "simulations": n_sim,
            "continuous_sign_test_power": float(reject_sign.mean()),
            "continuous_t_test_power": float(reject_t.mean()),
            "matched_sign_test_power": float(reject_match.mean()),
            "mean_estimated_slope": float(slope_mat.mean()),
            "continuous_heats": n_cont, "matched_heats": n_match,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=Path("results/q2/final_freeze_2026-11-28/final_curves.csv"))
    ap.add_argument("--output", type=Path, default=Path("results/q2/round4_identifiability_2026-11-28"))
    ap.add_argument("--simulations", type=int, default=N_SIM)
    ap.add_argument("--source-data-zip", type=Path, default=Path("data/raw/official/Data.zip"))
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw = args.input.read_bytes()
    curves = pd.read_csv(args.input)
    required = {
        "scenario", "heat", "time_s", "oxygen_cumulative_fit_m3",
        "oxygen_flow_analytic_m3_s", "CO2_flow_k2_m3_s", "CO_flow_k2_m3_s",
    }
    missing = required - set(curves.columns)
    if missing or len(curves) != 27178 or curves.heat.nunique() != 10:
        raise ValueError(f"Unexpected Q2.1 input: missing={sorted(missing)}, rows={len(curves)}")
    if curves[list(required)].isna().any().any():
        raise ValueError("Missing values in required Q2.1 columns")

    windows = estimate_windows(curves)
    windows.to_csv(args.output / "lag_windows.csv", index=False, float_format="%.10g")
    match_pair_rows, match_summary_rows, heat_slope_rows, pooled_rows, diag_rows = [], [], [], [], []
    for scenario in SCENARIOS:
        for gas in GASES:
            g = windows.query("scenario==@scenario and window_s==90 and metric=='pearson' and gas==@gas").copy()
            matched = round3_pairs(g, 0.15)
            for heat, h, pairs in matched:
                for p in pairs:
                    match_pair_rows.append({"scenario": scenario, "gas": gas, **p})
            for degree in [1, 2]:
                slopes = heat_slopes(g, degree)
                if len(slopes):
                    vals = slopes.slope_s_per_sd.to_numpy(float)
                    boot_rng = np.random.default_rng(SEED + degree + (0 if gas == "CO" else 100) + (0 if scenario == "H0" else 1000))
                    boot = boot_rng.choice(vals, size=(20000, len(vals)), replace=True).mean(axis=1)
                    p_t = float(ttest_1samp(vals, 0, alternative="less").pvalue) if len(vals) > 1 else 1.0
                    heat_slope_rows.append({
                        "scenario": scenario, "gas": gas, "stage_degree": degree,
                        "n_heats": len(vals), "mean_heat_slope_s_per_sd": float(vals.mean()),
                        "median_heat_slope_s_per_sd": float(np.median(vals)),
                        "heat_boot_ci_low": float(np.quantile(boot, 0.025)),
                        "heat_boot_ci_high": float(np.quantile(boot, 0.975)),
                        "exact_sign_p_beta_lt_0": exact_sign_p(vals), "one_sided_t_p_beta_lt_0": p_t,
                    })
                    for row in slopes.to_dict("records"):
                        diag_rows.append({"scenario": scenario, "gas": gas, "stage_degree": degree, **row})
            for degree in [1, 2]:
                for boundary in ["all", "interior"]:
                    result = pooled_cluster_regression(g, degree, boundary)
                    if result:
                        pooled_rows.append({"scenario": scenario, "gas": gas, **result})
            # Store matching summary using heat as the independent unit.
            heat_effects = []
            for heat, h, pairs in matched:
                diffs = [p["lag_diff_high_minus_low_s"] for p in pairs]
                heat_effects.append(float(np.median(diffs)))
            match_summary_rows.append({
                "scenario": scenario, "gas": gas,
                "n_pairs": sum(len(z[2]) for z in matched), "n_heats": len(matched),
                "median_heat_effect_s": float(np.median(heat_effects)) if heat_effects else np.nan,
                "exact_sign_p_high_flow_shorter": exact_sign_p(heat_effects),
            })

    pd.DataFrame(match_pair_rows).to_csv(args.output / "matched_pairs.csv", index=False, float_format="%.10g")
    pd.DataFrame(match_summary_rows).to_csv(args.output / "matching_summary.csv", index=False, float_format="%.10g")
    pd.DataFrame(heat_slope_rows).to_csv(args.output / "heat_slope_summary.csv", index=False, float_format="%.10g")
    pd.DataFrame(pooled_rows).to_csv(args.output / "pooled_cluster_regression.csv", index=False, float_format="%.10g")
    pd.DataFrame(diag_rows).to_csv(args.output / "heat_diagnostics.csv", index=False, float_format="%.10g")

    # Power analysis is restricted to the primary H0 / 90 s / Pearson estimand.
    power_linear = []
    power_quadratic = []
    for gas in GASES:
        g = windows.query("scenario=='H0' and window_s==90 and metric=='pearson' and gas==@gas").copy()
        matched = round3_pairs(g, 0.15)
        power_linear.append(power_simulation(g, 1, matched, args.simulations, SEED + (0 if gas == "CO" else 100)).assign(gas=gas))
        power_quadratic.append(power_simulation(g, 2, matched, args.simulations, SEED + (0 if gas == "CO" else 100) + 500).assign(gas=gas))
    pl = pd.concat(power_linear, ignore_index=True)
    pq = pd.concat(power_quadratic, ignore_index=True)
    pl.to_csv(args.output / "power_curve_linear_stage.csv", index=False, float_format="%.8g")
    pq.to_csv(args.output / "power_curve_quadratic_stage.csv", index=False, float_format="%.8g")

    mde_rows = []
    for label, frame in [("linear", pl), ("quadratic", pq)]:
        for gas, d in frame.groupby("gas"):
            for col in ["continuous_sign_test_power", "continuous_t_test_power", "matched_sign_test_power"]:
                hit = d[d[col] >= 0.80]
                mde = float(hit.effect_s_per_within_heat_sd.max()) if len(hit) else np.nan
                mde_rows.append({"stage_model": label, "gas": gas, "test": col,
                                 "mde_s_per_within_heat_sd": mde,
                                 "power_at_minus_16s": float(d.loc[d.effect_s_per_within_heat_sd == -16, col].iloc[0])})
    pd.DataFrame(mde_rows).to_csv(args.output / "mde_summary.csv", index=False)

    # Compact design diagnostics for the primary analysis.
    design_rows = []
    for gas in GASES:
        g = windows.query("scenario=='H0' and window_s==90 and metric=='pearson' and gas==@gas")
        for heat, h in g.groupby("heat"):
            p = h.oxygen_progress.to_numpy(float)
            q = h.flow_m3_s.to_numpy(float)
            X = stage_matrix(p, 1)
            xr, _ = residualize(q, X)
            design_rows.append({
                "gas": gas, "heat": heat, "n_windows": len(h),
                "flow_progress_corr": float(np.corrcoef(q, p)[0, 1]),
                "flow_resid_sd_after_linear_progress_m3_s": float(np.std(xr, ddof=1)),
                "boundary_peak_rate": float(h.boundary_peak.mean()),
            })
    pd.DataFrame(design_rows).to_csv(args.output / "design_diagnostics.csv", index=False, float_format="%.10g")

    source_zip_sha = None
    source_zip = args.source_data_zip
    if source_zip.exists():
        source_zip_sha = hashlib.sha256(source_zip.read_bytes()).hexdigest()
    qa = {
        "source_data_zip_sha256": source_zip_sha,
        "expected_source_data_zip_sha256": "cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727",
        "source_data_zip_hash_match": source_zip_sha == "cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727",
        "input_curves_sha256": hashlib.sha256(raw).hexdigest(),
        "repository_frozen_curves_reference_sha256": "9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49",
        "byte_hash_match_to_repository_reference": hashlib.sha256(raw).hexdigest() == "9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49",
        "input_rows": len(curves), "input_heats": int(curves.heat.nunique()),
        "input_scenarios": sorted(curves.scenario.unique().tolist()),
        "missing_cells_in_required_columns": int(curves[list(required)].isna().sum().sum()),
        "minimum_analytic_oxygen_flow_m3_s": float(curves.oxygen_flow_analytic_m3_s.min()),
        "oxygen_cumulative_monotone_all_heat_scenarios": bool(all(np.min(np.diff(g.oxygen_cumulative_fit_m3)) >= -1e-7 for _, g in curves.groupby(["scenario", "heat"]))),
        "h0_90s_pearson_median_lag_s": {gas: float(windows.query("scenario=='H0' and window_s==90 and metric=='pearson' and gas==@gas").lag_s.median()) for gas in GASES},
        "h0_round3_reproduced_match_median_effect_s": {r["gas"]: r["median_heat_effect_s"] for r in match_summary_rows if r["scenario"] == "H0"},
        "note": "Input curve byte hash differs from the repository-recorded frozen-curve hash; headline lag and round-three matched summaries are reproduced numerically, not byte-identically."
    }
    (args.output / "input_reconstruction_qa.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    manifest = {
        "experiment": "Q2.3 Round 4 — design-conditional power and identifiability audit",
        "input_curves_sha256": hashlib.sha256(raw).hexdigest(),
        "input_rows": len(curves), "window_rows": len(windows),
        "heats": HEATS, "scenarios": SCENARIOS, "windows_s": WINDOWS,
        "metrics": METRICS, "primary_estimand": "H0 / 90 s / Pearson",
        "flow_adjustment": "within-heat residual flow after linear or quadratic adjustment for cumulative oxygen progress; standardized within heat",
        "lag_grid_s": 2, "lag_search_range_s": "0 to window/6 (30 s for 90 s windows)",
        "power_effects_s_per_within_heat_sd": EFFECTS,
        "simulations_per_cell": args.simulations, "random_seed": SEED,
        "residual_resampling": "within-heat circular moving-block bootstrap, block length 2 windows",
        "inference": "heat-level exact sign test and one-sample t-test; matched contrasts use round-three pairs and heat-level exact sign test",
        "matched_design": "within-heat 25/75 flow tails, one-to-one Hungarian match on oxygen progress, 0.15 progress caliper",
        "source_data_zip_sha256": source_zip_sha,
        "limitations": [
            "Power simulation injects effects into estimated window-level lag outcomes, not raw physical signals.",
            "It is conditional on the empirical window design and residual distribution; it is not independent validation.",
            "The quadratic continuous analysis excludes heats with four or fewer windows from its per-heat slope test.",
            "Lag estimates at the search boundary are censored at the maximum searched lag.",
            "No causal effect is identified by this study."
        ],
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("Q2.3 ROUND 4 COMPLETE")
    print("Input rows:", len(curves), "lag windows:", len(windows))
    print("Output:", args.output)
    print("Source Data.zip SHA-256:", source_zip_sha)
    print("Input curves SHA-256:", manifest["input_curves_sha256"])
    print("Power simulation repetitions per cell:", args.simulations)


if __name__ == "__main__":
    main()
