#!/usr/bin/env python3
"""Prewhitened, windowed delay diagnostics for Q2.

Uses only the ten tracked process-curve CSVs. H0 keeps the official 0.5 s/row
clock. H1 is a sensitivity scenario that doubles the time axis and halves the
derived oxygen flow for heat 6C02057 only.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HEATS = [
    "6B03205", "6B03206", "6B03270", "6B03271", "6C02057",
    "6C02058", "6C02059", "6C02060", "6C02107", "6C02108",
]
METRICS = ("pearson", "spearman", "dcor")
SCENARIOS = ("H0", "H1")
WINDOWS = (90, 180)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xc, yc = x - x.mean(), y - y.mean()
    den = np.sqrt(np.dot(xc, xc) * np.dot(yc, yc))
    return float(np.dot(xc, yc) / den) if den > 0 else np.nan


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    ry = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    return pearson(rx, ry)


def distance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Biased, double-centred distance correlation used in Cattaldo et al. (2024)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    a = np.abs(x[:, None] - x[None, :])
    b = np.abs(y[:, None] - y[None, :])
    ac = a - a.mean(axis=1, keepdims=True) - a.mean(axis=0, keepdims=True) + a.mean()
    bc = b - b.mean(axis=1, keepdims=True) - b.mean(axis=0, keepdims=True) + b.mean()
    dcov2 = float(np.mean(ac * bc))
    dvar_x2 = float(np.mean(ac * ac))
    dvar_y2 = float(np.mean(bc * bc))
    den = np.sqrt(dvar_x2 * dvar_y2)
    return max(0.0, dcov2 / den) if den > 0 else 0.0


def detrend_linear(x: np.ndarray) -> np.ndarray:
    t = np.arange(len(x), dtype=float)
    design = np.column_stack((np.ones(len(x)), t))
    beta = np.linalg.lstsq(design, np.asarray(x, dtype=float), rcond=None)[0]
    return np.asarray(x, dtype=float) - design @ beta


def prewhiten(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Estimate AR(1) from detrended input and apply the same filter to output."""
    xd = detrend_linear(x)
    yd = detrend_linear(y)
    den = float(np.dot(xd[:-1], xd[:-1]))
    phi = float(np.dot(xd[:-1], xd[1:]) / den) if den > 0 else 0.0
    phi = float(np.clip(phi, -0.98, 0.98))
    xw = xd[1:] - phi * xd[:-1]
    yw = yd[1:] - phi * yd[:-1]
    return xw, yw, phi


def metric_value(x: np.ndarray, y: np.ndarray, metric: str) -> tuple[float, float, float]:
    r = pearson(x, y)
    if metric == "pearson":
        value = abs(r)
        signed = r
    elif metric == "spearman":
        signed = spearman(x, y)
        value = abs(signed)
    elif metric == "dcor":
        signed = distance_correlation(x, y)
        value = signed
    else:
        raise ValueError(metric)
    return float(value), float(signed), float(r)


def best_lag(
    x: np.ndarray, y: np.ndarray, max_lag_samples: int, metric: str
) -> tuple[int, float, float, float]:
    best = (-1, -np.inf, np.nan, np.nan)
    for k in range(max_lag_samples + 1):
        yy = y[k:k + len(x)]
        score, signed, r = metric_value(x, yy, metric)
        if np.isfinite(score) and score > best[1]:
            best = (k, score, signed, r)
    return best


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    m = len(p_values)
    order = np.argsort(p_values)
    q_sorted = np.empty(m, dtype=float)
    running = 1.0
    for j in range(m - 1, -1, -1):
        rank = j + 1
        running = min(running, p_values[order[j]] * m / rank)
        q_sorted[j] = min(1.0, running)
    q = np.empty(m, dtype=float)
    q[order] = q_sorted
    return q


def heat_weighted_fit(frame: pd.DataFrame) -> dict[str, float]:
    """Within-heat flow slope, with each heat receiving equal total weight."""
    groups = [g for _, g in frame.groupby("heat") if len(g)]
    if not groups:
        return {"slope": 0.0, "xbar": np.nan, "ybar": np.nan}
    heat_means = [
        (float(g.flow_m3_s.mean()), float(g.lag_s.mean())) for g in groups
    ]
    xbar = float(np.mean([v[0] for v in heat_means]))
    ybar = float(np.mean([v[1] for v in heat_means]))
    xx = 0.0
    xy = 0.0
    for g in groups:
        xm = float(g.flow_m3_s.mean())
        ym = float(g.lag_s.mean())
        xx += float(np.sum((g.flow_m3_s.to_numpy() - xm) ** 2)) / len(g)
        xy += float(np.sum(
            (g.flow_m3_s.to_numpy() - xm) * (g.lag_s.to_numpy() - ym)
        )) / len(g)
    return {"slope": xy / xx if xx > 0 else 0.0, "xbar": xbar, "ybar": ybar}


class LCG:
    """Small deterministic bootstrap RNG, matching the archived JS experiment."""
    def __init__(self, seed: int = 20261127):
        self.state = seed

    def index(self, n: int) -> int:
        self.state = (1664525 * self.state + 1013904223) & 0xFFFFFFFF
        return self.state % n


def bootstrap_slope(frame: pd.DataFrame, rng: LCG, n_boot: int = 5000) -> tuple[float, float]:
    per_heat = []
    for heat in HEATS:
        g = frame[frame.heat == heat]
        if len(g) < 2:
            per_heat.append((0.0, 0.0))
            continue
        xm = float(g.flow_m3_s.mean())
        ym = float(g.lag_s.mean())
        xx = float(np.sum((g.flow_m3_s.to_numpy() - xm) ** 2)) / len(g)
        xy = float(np.sum(
            (g.flow_m3_s.to_numpy() - xm) * (g.lag_s.to_numpy() - ym)
        )) / len(g)
        per_heat.append((xx, xy))
    estimates = []
    for _ in range(n_boot):
        xx = 0.0
        xy = 0.0
        for _ in HEATS:
            a, b = per_heat[rng.index(len(per_heat))]
            xx += a
            xy += b
        if xx > 0:
            estimates.append(xy / xx)
    estimates.sort()
    if not estimates:
        return np.nan, np.nan
    return estimates[int(0.025 * len(estimates))], estimates[int(0.975 * len(estimates))]


def make_tables(window_lags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries = []
    loho_rows = []
    flow_bins = []
    rng = LCG()

    for metric in METRICS:
        for scenario in SCENARIOS:
            for window_s in WINDOWS:
                for gas in ("CO2", "CO"):
                    g = window_lags[
                        (window_lags.metric == metric)
                        & (window_lags.scenario == scenario)
                        & (window_lags.window_s == window_s)
                        & (window_lags.gas == gas)
                    ].copy()
                    fit = heat_weighted_fit(g)
                    ci_low, ci_high = bootstrap_slope(g, rng)
                    loo_slopes = [
                        heat_weighted_fit(g[g.heat != heat])["slope"] for heat in HEATS
                    ]

                    test_mae = []
                    for test_heat in HEATS:
                        test = g[g.heat == test_heat]
                        train = g[g.heat != test_heat]
                        if test.empty or train.empty:
                            continue
                        train_fit = heat_weighted_fit(train)
                        fixed_pred = train_fit["ybar"]
                        flow_pred = (
                            train_fit["ybar"]
                            + train_fit["slope"] * (test.flow_m3_s.to_numpy() - train_fit["xbar"])
                        )
                        fixed_mae = float(np.mean(np.abs(test.lag_s.to_numpy() - fixed_pred)))
                        flow_mae = float(np.mean(np.abs(test.lag_s.to_numpy() - flow_pred)))
                        loho_rows.append({
                            "metric": metric, "scenario": scenario, "window_s": window_s,
                            "gas": gas, "test_heat": test_heat, "n_test": len(test),
                            "MAE_fixed_mean_s": fixed_mae, "MAE_flow_linear_s": flow_mae,
                            "gain_s": fixed_mae - flow_mae,
                        })
                        test_mae.append((fixed_mae, flow_mae))

                    flow = np.sort(g.flow_m3_s.to_numpy())
                    c1 = flow[(len(flow) - 1) // 3]
                    c2 = flow[(2 * (len(flow) - 1)) // 3]
                    groups = (
                        ("low", g[g.flow_m3_s <= c1]),
                        ("middle", g[(g.flow_m3_s > c1) & (g.flow_m3_s <= c2)]),
                        ("high", g[g.flow_m3_s > c2]),
                    )
                    for label, b in groups:
                        flow_bins.append({
                            "metric": metric, "scenario": scenario, "window_s": window_s,
                            "gas": gas, "flow_bin": label, "n_windows": len(b),
                            "median_flow_m3_s": float(b.flow_m3_s.median()),
                            "median_lag_s": float(b.lag_s.median()),
                        })

                    summaries.append({
                        "metric": metric, "scenario": scenario, "window_s": window_s,
                        "gas": gas, "n_windows": len(g), "n_heats": g.heat.nunique(),
                        "median_score": float(g.peak_score.median()),
                        "median_lag_s": float(g.lag_s.median()),
                        "BH_q_le_0_1_windows": int((g.bh_q <= 0.1).sum()),
                        "supported_heats": ";".join(sorted(g.loc[g.bh_q <= 0.1, "heat"].unique())),
                        "within_heat_slope_s_per_m3s": fit["slope"],
                        "cluster_bootstrap_95_low": ci_low,
                        "cluster_bootstrap_95_high": ci_high,
                        "LOO_slope_min": float(np.min(loo_slopes)),
                        "LOO_slope_max": float(np.max(loo_slopes)),
                        "LOHO_fixed_mean_MAE_s": float(np.mean([x[0] for x in test_mae])),
                        "LOHO_flow_linear_MAE_s": float(np.mean([x[1] for x in test_mae])),
                        "LOHO_gain_s": float(np.mean([x[0] - x[1] for x in test_mae])),
                        "LOHO_heats_improved": int(sum(x[1] < x[0] for x in test_mae)),
                    })
    return pd.DataFrame(summaries), pd.DataFrame(loho_rows), pd.DataFrame(flow_bins)


def run(output: Path) -> None:
    data_dir = ROOT / "data" / "processed" / "process_curves"
    frames = {}
    for heat in HEATS:
        path = data_dir / f"curves_{heat}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frames[heat] = pd.read_csv(path)

    rows = []
    for metric in METRICS:
        for scenario in SCENARIOS:
            for window_s in WINDOWS:
                input_n = window_s // 3
                max_lag_samples = window_s // 6
                for heat in HEATS:
                    frame = frames[heat]
                    is_fifth = heat == "6C02057"
                    time_scale = 2.0 if scenario == "H1" and is_fifth else 1.0
                    flow_scale = 0.5 if scenario == "H1" and is_fifth else 1.0
                    t = frame.time_s.to_numpy(dtype=float) * time_scale
                    q = frame.oxygen_flow_centered10.to_numpy(dtype=float) * flow_scale
                    for start_s in range(0, int(t[-1] - window_s) + 1, window_s):
                        grid = start_s + 2.0 * np.arange(window_s // 2)
                        x_raw = np.interp(grid[:input_n], t, q)
                        for gas in ("CO2", "CO"):
                            y_raw = np.interp(
                                grid, t, frame[f"{gas}_flow"].to_numpy(dtype=float)
                            )
                            xw, yw, phi = prewhiten(x_raw, y_raw)
                            peak_k, peak_score, association, r_at_peak = best_lag(
                                xw, yw, max_lag_samples, metric
                            )
                            null_scores = []
                            for shift in range(len(yw)):
                                if shift <= max_lag_samples or shift >= len(yw) - max_lag_samples:
                                    continue
                                shifted = np.roll(yw, -shift)
                                null_scores.append(
                                    best_lag(xw, shifted, max_lag_samples, metric)[1]
                                )
                            p_value = (
                                1 + sum(v >= peak_score for v in null_scores)
                            ) / (1 + len(null_scores))
                            rows.append({
                                "metric": metric, "scenario": scenario,
                                "window_s": window_s, "input_s": input_n * 2,
                                "max_lag_s": max_lag_samples * 2, "heat": heat,
                                "gas": gas, "start_s": start_s,
                                "flow_m3_s": float(np.mean(x_raw)),
                                "flow_sd_m3_s": float(np.std(x_raw)),
                                "ar1_phi": phi, "lag_s": peak_k * 2,
                                "peak_score": peak_score,
                                "association_at_peak": association,
                                "pearson_r_at_peak": r_at_peak,
                                "p_value": p_value, "n_null": len(null_scores),
                                "boundary_peak": bool(peak_k == max_lag_samples),
                            })

    window_lags = pd.DataFrame(rows)
    window_lags["bh_q"] = np.nan
    for _, idx in window_lags.groupby(
        ["metric", "scenario", "window_s", "gas"], sort=False
    ).groups.items():
        window_lags.loc[idx, "bh_q"] = bh_adjust(
            window_lags.loc[idx, "p_value"].to_numpy(dtype=float)
        )

    summary, loho, bins = make_tables(window_lags)
    output.mkdir(parents=True, exist_ok=True)
    window_lags.to_csv(output / "window_lags.csv", index=False, float_format="%.10g")
    summary.to_csv(output / "summary.csv", index=False, float_format="%.10g")
    loho.to_csv(output / "loho_scores.csv", index=False, float_format="%.10g")
    bins.to_csv(output / "flow_bins.csv", index=False, float_format="%.10g")
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", default=str(ROOT / "results" / "q2" / "prewhitened_delay" / "tables")
    )
    args = parser.parse_args()
    run(Path(args.output).resolve())


if __name__ == "__main__":
    main()
