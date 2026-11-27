#!/usr/bin/env python3
"""Compare native, regular, random, and bursty row-thinning mechanisms for 6C02057.

The target heat is never modified. Synthetic masks are applied only to the five
other 6C heats, which share the same native wavelength grid. Results are
exploratory and descriptive; five control heats do not support strong
population-level inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

TARGET_DEFAULT = "6C02057"
DT = 0.5


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def heat_id(path: Path) -> str | None:
    m = re.search(r"Heat No\.\s*(6[BC]\d+)\.xlsb$", path.name)
    return m.group(1) if m else None


def load_heat(path: Path) -> dict:
    frame = pd.read_excel(path, header=None, engine="pyxlsb")
    values = frame.iloc[1:].to_numpy(dtype=np.float64, copy=True)
    wl = frame.iloc[0, 1:2048].to_numpy(dtype=np.float64, copy=True)
    if values.shape[1] != 2051 or wl.size != 2047:
        raise ValueError(f"Unexpected workbook shape in {path.name}: {values.shape}")
    return {
        "heat": heat_id(path),
        "seq": values[:, 0].astype(np.int64),
        "X": values[:, 1:2048],
        "O": values[:, 2048],
        "CO2": values[:, 2049],
        "CO": values[:, 2050],
        "wl": wl,
        "path": path,
        "sha256": sha256(path),
        "missing_cells": int(np.isnan(values).sum()),
    }


def load_heat_cache(path: Path, source_path: Path | None = None) -> dict:
    z = np.load(path)
    required = {"seq", "X", "O", "CO2", "CO", "wl"}
    if not required.issubset(set(z.files)):
        raise ValueError(f"Cache {path.name} missing keys: {sorted(required - set(z.files))}")
    return {
        "heat": path.stem,
        "seq": z["seq"].astype(np.int64),
        "X": z["X"].astype(np.float64, copy=False),
        "O": z["O"].astype(np.float64, copy=False),
        "CO2": z["CO2"].astype(np.float64, copy=False),
        "CO": z["CO"].astype(np.float64, copy=False),
        "wl": z["wl"].astype(np.float64, copy=False),
        "path": source_path or path,
        "sha256": sha256(source_path) if source_path and source_path.exists() else sha256(path),
        "source_kind": "official_xlsb_with_npz_cache" if source_path and source_path.exists() else "npz_cache",
        "missing_cells": int(np.isnan(z["X"]).sum() + np.isnan(z["O"]).sum() + np.isnan(z["CO2"]).sum() + np.isnan(z["CO"]).sum()),
    }


def feature_series(d: dict) -> dict[str, np.ndarray]:
    X, wl = d["X"], d["wl"]
    row_sum = X.sum(axis=1)
    return {
        "CO2": d["CO2"],
        "CO": d["CO"],
        "spec_integral": np.trapezoid(X, wl, axis=1),
        "spec_centroid": (X @ wl) / np.maximum(row_sum, 1e-12),
        "spec_p95": np.quantile(X, 0.95, axis=1),
        "spec_max": X.max(axis=1),
        "spec_sat_fraction": np.mean(X == 65535, axis=1),
    }


def run_lengths(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x)
    if x.size == 0:
        return np.array([], dtype=int)
    starts = np.r_[0, np.flatnonzero(np.diff(x) != 0) + 1]
    ends = np.r_[starts[1:], len(x)]
    return ends - starts


def series_metrics(x: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    dx = np.diff(x)
    nz = np.abs(dx[dx != 0])
    runs = run_lengths(x)
    acf1 = float(np.corrcoef(x[:-1], x[1:])[0, 1]) if len(x) > 2 and np.std(x[:-1]) > 0 and np.std(x[1:]) > 0 else np.nan
    return {
        "n": int(len(x)),
        "transition_rate": float(np.mean(dx != 0)) if len(dx) else np.nan,
        "repeat_rate": float(np.mean(dx == 0)) if len(dx) else np.nan,
        "run_median": float(np.median(runs)) if len(runs) else np.nan,
        "run_p90": float(np.quantile(runs, 0.90)) if len(runs) else np.nan,
        "absdiff_median_all": float(np.median(np.abs(dx))) if len(dx) else np.nan,
        "absdiff_median_nonzero": float(np.median(nz)) if len(nz) else np.nan,
        "absdiff_p90_nonzero": float(np.quantile(nz, 0.90)) if len(nz) else np.nan,
        "acf1": acf1,
        "mean": float(np.mean(x)),
        "sd": float(np.std(x)),
    }


def oxygen_metrics(O: np.ndarray, quantum: float) -> dict[str, float]:
    O = np.asarray(O, dtype=float)
    delta = np.diff(O)
    positive = delta[delta > 0]
    update_rows = np.flatnonzero(delta > 0) + 1
    gaps = np.diff(np.r_[0, update_rows])
    if positive.size == 0:
        return {k: np.nan for k in ["flat_rate", "positive_update_rate", "gap1_share", "gap2_share", "gap3plus_share", "increment_median", "increment_mean", "increment_p90", "multi_increment_share"]}
    return {
        "flat_rate": float(np.mean(delta == 0)),
        "positive_update_rate": float(np.mean(delta > 0)),
        "gap1_share": float(np.mean(gaps == 1)),
        "gap2_share": float(np.mean(gaps == 2)),
        "gap3plus_share": float(np.mean(gaps >= 3)),
        "increment_median": float(np.median(positive)),
        "increment_mean": float(np.mean(positive)),
        "increment_p90": float(np.quantile(positive, 0.90)),
        "multi_increment_share": float(np.mean(positive > 1.5 * quantum)),
    }


def mask_regular(n: int, m: int, rng: np.random.Generator) -> np.ndarray:
    # Systematic near-equally-spaced retention with a randomized phase.
    phase = rng.random()
    return np.floor((phase + np.arange(m)) * n / m).astype(int)


def mask_markov(n: int, m: int, mean_cycle: float, rng: np.random.Generator) -> np.ndarray:
    """Alternating geometric keep/drop runs, adjusted to exactly m retained rows."""
    p_keep = m / n
    mean_keep = max(1.01, mean_cycle * p_keep)
    mean_drop = max(1.01, mean_cycle * (1.0 - p_keep))
    p_kd = min(1.0, 1.0 / mean_keep)
    p_dk = min(1.0, 1.0 / mean_drop)
    state = bool(rng.random() < p_keep)
    keep = np.zeros(n, dtype=bool)
    for i in range(n):
        keep[i] = state
        if state:
            if rng.random() < p_kd:
                state = False
        else:
            if rng.random() < p_dk:
                state = True
    diff = int(keep.sum()) - m
    if diff > 0:
        candidates = np.flatnonzero(keep)
        keep[rng.choice(candidates, size=diff, replace=False)] = False
    elif diff < 0:
        candidates = np.flatnonzero(~keep)
        keep[rng.choice(candidates, size=-diff, replace=False)] = True
    return np.flatnonzero(keep)


def make_mask(n: int, m: int, scenario: str, rng: np.random.Generator) -> np.ndarray:
    if scenario == "regular":
        return mask_regular(n, m, rng)
    if scenario == "random":
        return np.sort(rng.choice(n, size=m, replace=False))
    if scenario.startswith("bursty"):
        mean_cycle = float(scenario.replace("bursty", ""))
        return mask_markov(n, m, mean_cycle, rng)
    raise ValueError(scenario)


def feature_metrics_for_sample(d: dict, idx: np.ndarray, quantum: float) -> list[dict]:
    rows = []
    for sensor, x in d["features"].items():
        for metric, value in series_metrics(x[idx]).items():
            if metric == "n":
                continue
            rows.append({"sensor": sensor, "feature": metric, "value": value})
    for metric, value in oxygen_metrics(d["O"][idx], quantum).items():
        rows.append({"sensor": "oxygen_counter", "feature": metric, "value": value})
    return rows


def cross_corr(x: np.ndarray, y: np.ndarray, max_lag: int = 10) -> tuple[float, int, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    vals = []
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            a, b = x[-lag:], y[: n + lag]
        elif lag > 0:
            a, b = x[: n - lag], y[lag:]
        else:
            a, b = x, y
        if len(a) > 5 and np.std(a) > 0 and np.std(b) > 0:
            vals.append((float(np.corrcoef(a, b)[0, 1]), lag))
    if not vals:
        return np.nan, 0, np.nan
    corr0 = next((v for v, lag in vals if lag == 0), np.nan)
    best = max(vals, key=lambda z: abs(z[0]))
    return corr0, int(best[1]), float(best[0])


def cross_sensor_features(d: dict, idx_o: np.ndarray, idx_g: np.ndarray, idx_s: np.ndarray, threshold: float) -> dict[str, float]:
    O = d["O"][idx_o]
    CO2 = d["CO2"][idx_g]
    CO = d["CO"][idx_g]
    spec = d["features"]["spec_integral"][idx_s]
    o_flag = (np.diff(O) > 0).astype(float)
    c2_flag = (np.diff(CO2) != 0).astype(float)
    co_flag = (np.diff(CO) != 0).astype(float)
    s_flag = (np.abs(np.diff(spec)) > threshold).astype(float)
    out = {}
    for name, gas_flag in [("co2", c2_flag), ("co", co_flag), ("spectrum", s_flag)]:
        c0, lag, cmax = cross_corr(o_flag, gas_flag, 10)
        out[f"o2_{name}_corr0"] = c0
        out[f"o2_{name}_best_lag"] = lag
        out[f"o2_{name}_best_corr"] = cmax
    out["co2_co_corr0"] = cross_corr(c2_flag, co_flag, 0)[0]
    return out


def spectral_shape_metrics(X: np.ndarray, idx: np.ndarray) -> dict[str, float]:
    A, B = X[idx[:-1]], X[idx[1:]]
    norm_a = np.linalg.norm(A, axis=1)
    norm_b = np.linalg.norm(B, axis=1)
    cosine = np.sum(A * B, axis=1) / np.maximum(norm_a * norm_b, 1e-12)
    Ac = A - A.mean(axis=1, keepdims=True)
    Bc = B - B.mean(axis=1, keepdims=True)
    pearson = np.sum(Ac * Bc, axis=1) / np.maximum(np.linalg.norm(Ac, axis=1) * np.linalg.norm(Bc, axis=1), 1e-12)
    relrmse = np.sqrt(np.mean((B - A) ** 2, axis=1)) / np.maximum(np.sqrt(np.mean(A ** 2, axis=1)), 1e-12)
    return {
        "cosine_median": float(np.median(cosine)), "cosine_p10": float(np.quantile(cosine, .1)),
        "pearson_median": float(np.median(pearson)), "pearson_p10": float(np.quantile(pearson, .1)),
        "relrmse_median": float(np.median(relrmse)), "relrmse_p90": float(np.quantile(relrmse, .9)),
    }


def empirical_percentile(values: np.ndarray, target: float) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values) or not np.isfinite(target):
        return np.nan
    less = np.sum(values < target)
    equal = np.sum(values == target)
    return float((less + 0.5 * equal) / len(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, help="Directory containing official XLSB workbooks")
    parser.add_argument("--cache-dir", type=Path, help="Optional directory of pre-extracted per-heat NPZ files")
    parser.add_argument("--out", type=Path, default=Path("results/q2/missingness_mechanism"))
    parser.add_argument("--target", default=TARGET_DEFAULT)
    parser.add_argument("--n-rep", type=int, default=250)
    parser.add_argument("--cross-rep", type=int, default=200)
    parser.add_argument("--shape-rep", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261128)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    if args.cache_dir:
        cache_paths = sorted(p for p in args.cache_dir.glob("*.npz") if p.stem.startswith("6C"))
        data = {}
        for cp in cache_paths:
            source = None
            if args.data_dir:
                source = next((p for p in args.data_dir.glob("*.xlsb") if heat_id(p) == cp.stem), None)
            data[cp.stem] = load_heat_cache(cp, source)
    else:
        if not args.data_dir:
            parser.error("Provide --data-dir or --cache-dir")
        paths = sorted(p for p in args.data_dir.glob("*.xlsb") if heat_id(p) and heat_id(p).startswith("6C"))
        data = {heat_id(p): load_heat(p) for p in paths}
    if args.target not in data or len(data) != 6:
        raise ValueError(f"Expected the six 6C heats including {args.target}; found {sorted(data)}")
    target = data[args.target]
    controls = {h: d for h, d in data.items() if h != args.target}
    for d in data.values():
        d["features"] = feature_series(d)
    grids = [d["wl"] for d in data.values()]
    if not all(np.array_equal(grids[0], w) for w in grids[1:]):
        raise ValueError("6C wavelength grids are not identical; harmonization is required")

    # The 6C positive-increment median is used only as a descriptive counter quantum proxy.
    positive_all = np.concatenate([np.diff(d["O"])[np.diff(d["O"]) > 0] for d in data.values()])
    quantum = float(np.median(positive_all))
    control_flows = {h: (d["O"][-1] - d["O"][0]) / ((len(d["O"]) - 1) * DT) for h, d in controls.items()}
    reference_flow = float(np.median(list(control_flows.values())))
    target_net = float(target["O"][-1] - target["O"][0])
    estimated_full_duration = target_net / reference_flow
    estimated_full_rows = estimated_full_duration / DT + 1.0
    p_hat = float(len(target["O"]) / estimated_full_rows)
    target_n = len(target["O"])

    audit_rows = []
    for h, d in data.items():
        audit_rows.append({
            "heat": h, "rows": len(d["seq"]),
            "sequence_contiguous": bool(np.array_equal(d["seq"], np.arange(1, len(d["seq"]) + 1))),
            "missing_cells": d["missing_cells"], "sha256": d["sha256"],
            "source_kind": d.get("source_kind", "official_xlsb"),
            "oxygen_net_m3": float(d["O"][-1] - d["O"][0]),
            "mean_flow_at_0p5_m3_s": float((d["O"][-1] - d["O"][0]) / ((len(d["O"]) - 1) * DT)),
        })
    pd.DataFrame(audit_rows).to_csv(args.out / "data_audit.csv", index=False)

    # Target feature vector under its observed row order.
    target_rows = []
    for row in feature_metrics_for_sample(target, np.arange(target_n), quantum):
        target_rows.append({"heat": args.target, **row})
    pd.DataFrame(target_rows).to_csv(args.out / "target_metrics.csv", index=False)

    rng = np.random.default_rng(args.seed)
    scenarios = ["regular", "random", "bursty4", "bursty8", "bursty16", "bursty32"]
    simulation_rows = []
    shape_rows = []
    cross_rows = []

    shape_rows.append({"heat": args.target, "scenario": "target_native", "rep": 0,
                       **spectral_shape_metrics(target["X"], np.arange(target_n))})

    for h, d in controls.items():
        n = len(d["O"])
        shape_rows.append({"heat": h, "scenario": "native", "rep": 0,
                           **spectral_shape_metrics(d["X"], np.arange(n))})
        m_fixed = min(target_n, n - 1)
        m_rate = max(2, min(n - 1, int(round(n * p_hat))))
        native_idx = np.arange(n)
        for row in feature_metrics_for_sample(d, native_idx, quantum):
            simulation_rows.append({"heat": h, "scenario": "native", "rep": 0,
                                    "retention_rate": 1.0, **row})

        # Exact target-row-count experiment: all scenarios retain the same number of rows.
        for scenario in scenarios:
            for rep in range(args.n_rep):
                idx = make_mask(n, m_fixed, scenario, rng)
                for row in feature_metrics_for_sample(d, idx, quantum):
                    simulation_rows.append({"heat": h, "scenario": f"{scenario}_fixed_n", "rep": rep,
                                            "retention_rate": len(idx) / n, **row})

        # Retention-rate sensitivity around the oxygen-anchored estimate.
        for scenario in ["random", "bursty8", "bursty16", "bursty32"]:
            for rep in range(args.n_rep):
                idx = make_mask(n, m_rate, scenario, rng)
                for row in feature_metrics_for_sample(d, idx, quantum):
                    simulation_rows.append({"heat": h, "scenario": f"{scenario}_p_hat", "rep": rep,
                                            "retention_rate": len(idx) / n, **row})

        # Full-spectrum adjacent-frame similarity, sampled at lower replication count.
        X = d["X"]
        for scenario in ["regular", "random", "bursty8", "bursty32"]:
            for rep in range(min(args.shape_rep, args.n_rep)):
                idx = make_mask(n, m_fixed, scenario, rng)
                shape_rows.append({"heat": h, "scenario": f"{scenario}_fixed_n", "rep": rep,
                                   **spectral_shape_metrics(X, idx)})

        # Cross-sensor alignment: global row deletion versus independently compressed streams.
        spec_integral = d["features"]["spec_integral"]
        spec_threshold = float(np.quantile(np.abs(np.diff(spec_integral)), .75))
        for scenario in ["global", "independent"]:
            for rep in range(args.cross_rep):
                if scenario == "global":
                    idx = make_mask(n, m_fixed, "random", rng)
                    cross = cross_sensor_features(d, idx, idx, idx, spec_threshold)
                else:
                    idx_o = make_mask(n, m_fixed, "random", rng)
                    idx_g = make_mask(n, m_fixed, "random", rng)
                    idx_s = make_mask(n, m_fixed, "random", rng)
                    cross = cross_sensor_features(d, idx_o, idx_g, idx_s, spec_threshold)
                cross_rows.append({"heat": h, "scenario": scenario, "rep": rep, **cross})

    # Target cross-sensor synchrony under its observed row pairing.
    target_spec = target["features"]["spec_integral"]
    target_threshold = float(np.quantile(np.abs(np.diff(target_spec)), .75))
    target_idx = np.arange(target_n)
    target_cross = cross_sensor_features(target, target_idx, target_idx, target_idx, target_threshold)
    cross_rows.append({"heat": args.target, "scenario": "target", "rep": 0, **target_cross})

    sim = pd.DataFrame(simulation_rows)
    sim.to_csv(args.out / "simulation_metrics.csv", index=False)
    shape_df = pd.DataFrame(shape_rows)
    shape_df.to_csv(args.out / "spectral_shape_simulations.csv", index=False)
    cross = pd.DataFrame(cross_rows)
    cross.to_csv(args.out / "cross_sensor_simulations.csv", index=False)

    # Compare the target with each control heat's simulated distribution. Monte Carlo
    # replicates are summarized within heat first; they are not treated as extra heats.
    tdf = pd.DataFrame(target_rows)
    match_rows = []
    control_match_rows = []
    for (sensor, feature), tg in tdf.groupby(["sensor", "feature"]):
        target_value = float(tg.iloc[0]["value"])
        subset = sim[(sim.sensor == sensor) & (sim.feature == feature)]
        for scenario, ss in subset.groupby("scenario"):
            heat_scores, robust_scores, in80 = [], [], []
            for h, hs in ss.groupby("heat"):
                vals = hs.value.dropna().to_numpy(float)
                if not len(vals):
                    continue
                pct = empirical_percentile(vals, target_value)
                med = float(np.median(vals))
                q10, q25, q75, q90 = np.quantile(vals, [.1, .25, .75, .9])
                scale = max(float(q75 - q25), float((q90 - q10) / 2.563),
                            .05 * max(abs(med), abs(target_value), 1e-6))
                robust_z = abs(target_value - med) / scale
                inside = float(q10 <= target_value <= q90)
                heat_scores.append(abs(2 * pct - 1))
                robust_scores.append(robust_z)
                in80.append(inside)
                control_match_rows.append({"sensor": sensor, "feature": feature, "scenario": scenario,
                                           "control_heat": h, "target_value": target_value,
                                           "simulation_median": med, "simulation_p10": float(q10),
                                           "simulation_p90": float(q90), "robust_scale": scale,
                                           "robust_distance": robust_z, "target_inside_p10_p90": inside,
                                           "target_percentile": pct})
            if heat_scores:
                match_rows.append({"sensor": sensor, "feature": feature, "scenario": scenario,
                                   "target_value": target_value,
                                   "median_heat_percentile_distance": float(np.median(heat_scores)),
                                   "median_heat_robust_distance": float(np.median(robust_scores)),
                                   "controls_target_inside_p10_p90_share": float(np.mean(in80)),
                                   "n_control_heats": len(heat_scores)})
    match = pd.DataFrame(match_rows)
    match.to_csv(args.out / "scenario_match_by_feature.csv", index=False)
    pd.DataFrame(control_match_rows).to_csv(args.out / "scenario_match_by_control.csv", index=False)

    # Full-spectrum shape similarity is kept separate from scalar spectral features.
    target_shape = shape_df[shape_df.scenario == "target_native"].iloc[0]
    shape_metrics = ["cosine_median", "cosine_p10", "pearson_median", "pearson_p10", "relrmse_median", "relrmse_p90"]
    shape_match_rows = []
    for scenario, ss in shape_df[shape_df.heat != args.target].groupby("scenario"):
        for metric in shape_metrics:
            per_heat = ss.groupby("heat")[metric].median().dropna()
            values = per_heat.to_numpy(float)
            target_value = float(target_shape[metric])
            if not len(values):
                continue
            med = float(np.median(values))
            scale = max(float(np.quantile(values, .75) - np.quantile(values, .25)),
                        float((np.quantile(values, .9) - np.quantile(values, .1)) / 2.563),
                        .05 * max(abs(med), abs(target_value), 1e-6))
            shape_match_rows.append({"scenario": scenario, "metric": metric, "target_value": target_value,
                                     "control_median": med,
                                     "robust_distance": abs(target_value - med) / scale,
                                     "controls_target_inside_p10_p90": bool(np.quantile(values, .1) <= target_value <= np.quantile(values, .9)),
                                     "n_control_heats": int(len(values))})
    pd.DataFrame(shape_match_rows).to_csv(args.out / "spectral_shape_match.csv", index=False)

    family_map = {
        "CO2": ["CO2"], "CO": ["CO"], "oxygen_counter": ["oxygen_counter"],
        "spectrum_scalar": ["spec_integral", "spec_centroid", "spec_p95", "spec_max", "spec_sat_fraction"],
    }
    summary_rows = []
    for family, sensors in family_map.items():
        for scenario, ss in match[match.sensor.isin(sensors)].groupby("scenario"):
            summary_rows.append({"family": family, "scenario": scenario,
                                 "median_percentile_distance": float(ss.median_heat_percentile_distance.median()),
                                 "median_robust_distance": float(ss.median_heat_robust_distance.median()),
                                 "mean_robust_distance": float(ss.median_heat_robust_distance.mean()),
                                 "mean_controls_inside_p10_p90": float(ss.controls_target_inside_p10_p90_share.mean()),
                                 "n_features": int(len(ss))})
    pd.DataFrame(summary_rows).to_csv(args.out / "scenario_match_summary.csv", index=False)

    cross_summary = []
    for (scenario, metric), ss in cross[cross.heat != args.target].melt(
        id_vars=["heat", "scenario", "rep"], var_name="metric", value_name="value"
    ).groupby(["scenario", "metric"]):
        per_heat = ss.groupby("heat").value.median().dropna()
        cross_summary.append({"scenario": scenario, "metric": metric,
                              "median_across_control_heats": float(per_heat.median()) if len(per_heat) else np.nan,
                              "control_heat_min": float(per_heat.min()) if len(per_heat) else np.nan,
                              "control_heat_max": float(per_heat.max()) if len(per_heat) else np.nan,
                              "n_control_heats": int(len(per_heat))})
    for metric, value in target_cross.items():
        cross_summary.append({"scenario": "target", "metric": metric,
                              "median_across_control_heats": float(value),
                              "control_heat_min": np.nan, "control_heat_max": np.nan,
                              "n_control_heats": 0})
    pd.DataFrame(cross_summary).to_csv(args.out / "cross_sensor_summary.csv", index=False)

    run_summary = {
        "target": args.target,
        "target_rows": target_n,
        "target_oxygen_net_m3": target_net,
        "control_mean_flows_m3_s": control_flows,
        "reference_flow_median_m3_s": reference_flow,
        "estimated_full_duration_s_if_target_flow_equals_reference_median": estimated_full_duration,
        "estimated_full_rows_at_0p5s_if_target_flow_equals_reference_median": estimated_full_rows,
        "implied_retention_rate": p_hat,
        "counter_increment_reference_quantum_proxy_m3": quantum,
        "controls": sorted(controls),
        "n_rep": args.n_rep,
        "cross_rep": args.cross_rep,
        "shape_rep": args.shape_rep,
        "seed": args.seed,
        "interpretation": "Exploratory mechanism comparison only; not proof of a particular acquisition failure mode.",
    }
    (args.out / "summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
    print(json.dumps(run_summary, indent=2))
    print("Wrote results to", args.out.resolve())


if __name__ == "__main__":
    main()
