#!/usr/bin/env python3
"""Reproduce Q3 physical-mechanism diagnostics from the official Data.zip.

This script reproduces the tail, retrospective stage, stage-lag, and spectral
transition diagnostics. It intentionally does not reproduce the pilot Ridge
ablation; that runner is tracked separately as pending work in the report.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyxlsb import open_workbook
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import detrend
from scipy.stats import theilslopes, wilcoxon

DT = 0.5


def load_data_zip(path: Path):
    heats = {}
    with zipfile.ZipFile(path) as archive:
        names = sorted(n for n in archive.namelist() if n.lower().endswith(".xlsb"))
        for name in names:
            with open_workbook(io.BytesIO(archive.read(name))) as wb:
                sheet_name = wb.sheets[0]
                with wb.get_sheet(sheet_name) as sheet:
                    rows = sheet.rows()
                    header = [cell.v for cell in next(rows)]
                    values = np.asarray([[cell.v for cell in row] for row in rows], dtype=float)
            if values.shape[1] != 2051:
                raise ValueError(f"Unexpected column count in {name}: {values.shape}")
            heat = str(sheet_name)
            heats[heat] = {
                "wavelength": np.asarray(header[1:2048], dtype=float),
                "spectrum": values[:, 1:2048],
                "oxygen": values[:, 2048],
                "gas": values[:, 2049:2051],  # CO2, CO; m3/s
                "file": name,
            }
    if len(heats) != 10:
        raise ValueError(f"Expected ten heats, got {len(heats)}")
    return heats


def smooth_gas(gas):
    # 11 rows ~= 5.5 s on the official 0.5 s grid; smooth time only, not gas channels.
    return median_filter(np.asarray(gas, dtype=float), size=(11, 1), mode="nearest")


def terminal_cutoff(data):
    gas = smooth_gas(data["gas"])
    n = len(gas)
    t = np.arange(n) * DT
    peak = np.max(gas, axis=0)
    low = np.all(gas < peak[None, :] * 0.05, axis=1)
    run = int(10 / DT)
    start_search = max(int(120 / DT), int(0.5 * n))
    for i in range(start_search, n - run):
        if low[i:i + run].all():
            end = min(i + run, n - 1)
            if data["oxygen"][end] - data["oxygen"][i] > 1.0:
                return float(t[i])
    return None


def robust_slope_per_min(t, y):
    if len(t) < 5:
        return np.nan
    return float(theilslopes(y, t).slope * 60.0)


def tail_diagnostics(heats, output_dir):
    rows = []
    sensitivity = []
    for heat, data in heats.items():
        gas = smooth_gas(data["gas"])
        t = np.arange(len(gas)) * DT
        co2, co = gas.T
        total = co2 + co
        share = co2 / np.maximum(total, 1e-6)
        cut = terminal_cutoff(data)
        end = cut if cut is not None else float(t[-1])
        for window in (30, 45, 60, 90):
            mask = (t >= end - window) & (t <= end)
            slopes = [robust_slope_per_min(t[mask], y[mask]) for y in (co2, co, total, share)]
            sensitivity.append([heat, window, *slopes])
        mask = (t >= end - 60) & (t <= end)
        slopes = [robust_slope_per_min(t[mask], y[mask]) for y in (co2, co, total, share)]

        def trailing_median(arr, at):
            m = (t >= at - 10) & (t <= at)
            return float(np.median(arr[m])) if m.any() else np.nan

        start = end - 60
        start_vals = [trailing_median(y, start) for y in (co2, co, total, share)]
        end_vals = [trailing_median(y, end) for y in (co2, co, total, share)]
        rows.append([
            heat, end, cut is not None, *slopes,
            start_vals[0], end_vals[0], start_vals[1], end_vals[1],
            start_vals[2], end_vals[2], start_vals[3], end_vals[3],
            end_vals[3] - start_vals[3],
        ])
    cols = ["heat", "tail_end_s", "terminal_collapse", "slope_CO2_per_min", "slope_CO_per_min",
            "slope_sum_per_min", "slope_share_per_min", "CO2_start", "CO2_end", "CO_start", "CO_end",
            "sum_start", "sum_end", "share_start", "share_end", "share_delta"]
    pd.DataFrame(rows, columns=cols).to_csv(output_dir / "tail_dynamics.csv", index=False)
    pd.DataFrame(sensitivity, columns=["heat", "window_s", "slope_CO2_per_min", "slope_CO_per_min",
                                      "slope_sum_per_min", "slope_share_per_min"]).to_csv(
        output_dir / "tail_window_sensitivity.csv", index=False)


def segment_cost_matrix(y, t, min_points):
    n, p = y.shape
    st = np.r_[0.0, np.cumsum(t)]
    st2 = np.r_[0.0, np.cumsum(t * t)]
    sy = np.vstack([np.zeros(p), np.cumsum(y, axis=0)])
    sty = np.vstack([np.zeros(p), np.cumsum(t[:, None] * y, axis=0)])
    sy2 = np.vstack([np.zeros(p), np.cumsum(y * y, axis=0)])
    cost = np.full((n + 1, n + 1), np.inf)
    for i in range(n):
        for j in range(i + min_points, n + 1):
            count = j - i
            sum_t = st[j] - st[i]
            sum_t2 = st2[j] - st2[i]
            det = count * sum_t2 - sum_t * sum_t
            sum_y = sy[j] - sy[i]
            sum_ty = sty[j] - sty[i]
            sum_y2 = sy2[j] - sy2[i]
            if det <= 1e-12:
                cost[i, j] = float(np.sum(sum_y2))
                continue
            b0 = (sum_t2 * sum_y - sum_t * sum_ty) / det
            b1 = (count * sum_ty - sum_t * sum_y) / det
            cost[i, j] = float(np.sum(sum_y2 - b0 * sum_y - b1 * sum_ty))
    return cost


def piecewise_linear_dp(y, t, segments=4, min_points=6):
    n = len(t)
    cost = segment_cost_matrix(y, t, min_points)
    dp = np.full((segments + 1, n + 1), np.inf)
    prev = np.full((segments + 1, n + 1), -1, dtype=int)
    dp[0, 0] = 0.0
    for k in range(1, segments + 1):
        for j in range(k * min_points, n + 1):
            starts = np.arange((k - 1) * min_points, j - min_points + 1)
            vals = dp[k - 1, starts] + cost[starts, j]
            if len(vals) and np.isfinite(vals).any():
                q = int(np.argmin(vals))
                dp[k, j] = vals[q]
                prev[k, j] = starts[q]
    bounds = [n]
    j = n
    for k in range(segments, 0, -1):
        i = prev[k, j]
        if i < 0:
            raise RuntimeError("No valid piecewise segmentation")
        bounds.append(int(i))
        j = int(i)
    return bounds[::-1]


def stage_boundaries(heat, data):
    gas = smooth_gas(data["gas"])
    t = np.arange(len(gas)) * DT
    cut = terminal_cutoff(data)
    n = len(gas) if cut is None else max(1, int(cut / DT))
    gas = gas[:n]
    t = t[:n]
    co2, co = gas.T
    total = co2 + co
    share = co2 / np.maximum(total, 0.1)
    share = np.clip(share, 0.01, 0.99)
    y = np.column_stack([np.log(total + 0.2), np.log(share / (1.0 - share))])
    y = gaussian_filter1d(y, sigma=8 / DT, axis=0, mode="nearest")
    idx = np.arange(0, len(y), int(5 / DT))
    y = y[idx]
    t5 = t[idx]
    y = (y - y.mean(axis=0)) / np.maximum(y.std(axis=0), 1e-9)
    bounds = piecewise_linear_dp(y, t5, segments=4, min_points=6)
    return [float(t5[i]) for i in bounds[1:-1]], n


def stage_diagnostics(heats, output_dir):
    bound_map = {}
    bound_rows = []
    slope_rows = []
    for heat, data in heats.items():
        bounds, n = stage_boundaries(heat, data)
        bound_map[heat] = bounds
        bound_rows.append([heat, *bounds, (n - 1) * DT])
        gas = smooth_gas(data["gas"])[:n]
        t = np.arange(n) * DT
        co2, co = gas.T
        total = co2 + co
        share = np.clip(co2 / np.maximum(total, 0.1), 0.01, 0.99)
        values = np.column_stack([
            np.log(co2 + 0.2), np.log(co + 0.2), np.log(total + 0.2),
            np.log(share / (1.0 - share)),
        ])
        values = gaussian_filter1d(values, sigma=8 / DT, axis=0, mode="nearest")
        indices = [0] + [int(round(x / DT)) for x in bounds] + [n]
        for stage, (i, j) in enumerate(zip(indices[:-1], indices[1:]), start=1):
            slopes = [robust_slope_per_min(t[i:j], values[i:j, k]) for k in range(4)]
            slope_rows.append([heat, stage, t[i], t[j - 1], *slopes])
    pd.DataFrame(bound_rows, columns=["heat", "b1_s", "b2_s", "b3_s", "analysis_end_s"]).to_csv(
        output_dir / "stage_boundaries.csv", index=False)
    slope_df = pd.DataFrame(slope_rows, columns=["heat", "stage", "start_s", "end_s", "logco2_slope",
                                                 "logco_slope", "logtotal_slope", "logitshare_slope"])
    slope_df.to_csv(output_dir / "stage_slopes.csv", index=False)
    summary = slope_df.groupby("stage").agg(
        n=("heat", "nunique"),
        co2_positive=("logco2_slope", lambda x: int((x > 0).sum())),
        median_logco2_slope=("logco2_slope", "median"),
        co_positive=("logco_slope", lambda x: int((x > 0).sum())),
        median_logco_slope=("logco_slope", "median"),
        total_positive=("logtotal_slope", lambda x: int((x > 0).sum())),
        median_logtotal_slope=("logtotal_slope", "median"),
        share_positive=("logitshare_slope", lambda x: int((x > 0).sum())),
        median_logitshare_slope=("logitshare_slope", "median"),
    ).reset_index()
    summary.to_csv(output_dir / "stage_summary.csv", index=False)
    return bound_map


def spectral_brightness(data):
    wl, x = data["wavelength"], data["spectrum"]
    edges = np.linspace(345, 1015, 33)
    bands = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (wl >= lo) & (wl < hi)
        bands.append(np.log1p(np.mean(x[:, m], axis=1)))
    return np.mean(np.column_stack(bands), axis=1)


def stage_lag_diagnostics(heats, bounds, output_dir):
    rows = []
    lags = np.arange(0, 61, 5)
    for heat, data in heats.items():
        brightness = spectral_brightness(data)
        gas = smooth_gas(data["gas"])
        n = len(gas)
        stage_times = bounds[heat]
        indices = [0] + [int(round(x / DT)) for x in stage_times] + [n]
        for stage, (i0, i1) in enumerate(zip(indices[:-1], indices[1:]), start=1):
            for gas_idx, gas_name in enumerate(("CO2", "CO")):
                profile = []
                for lag in lags:
                    d = int(lag / DT)
                    det = np.arange(max(i0 + d + 10, 11), min(i1, n), 2)
                    src = det - d
                    if len(det) < 20:
                        profile.append(0.0)
                        continue
                    dx = detrend(brightness[src] - brightness[src - 10])
                    dy = detrend(gas[det, gas_idx] - gas[det - 10, gas_idx])
                    corr = np.corrcoef(dx, dy)[0, 1] if dx.std() > 1e-12 and dy.std() > 1e-12 else 0.0
                    profile.append(abs(float(corr)) if np.isfinite(corr) else 0.0)
                best = int(lags[int(np.argmax(profile))])
                rows.append([heat, gas_name, stage, i0 * DT, i1 * DT, best, max(profile)])
    detail = pd.DataFrame(rows, columns=["heat", "gas", "stage", "start_s", "end_s", "best_lag_s", "abs_corr"])
    detail.to_csv(output_dir / "stage_lags_by_heat.csv", index=False)
    summary = detail.groupby(["gas", "stage"]).agg(
        n=("heat", "nunique"), median_lag_s=("best_lag_s", "median"),
        mean_abs_corr=("abs_corr", "mean"), median_abs_corr=("abs_corr", "median"),
    ).reset_index()
    summary.to_csv(output_dir / "stage_lag_summary.csv", index=False)


def band_mean_log(x, wl, lo, hi):
    mask = (wl >= lo) & (wl < hi)
    xx = x[:, mask].astype(float).copy()
    xx[xx >= 65535] = np.nan
    return np.nanmean(np.log1p(xx), axis=1)


def local_contrast(x, wl, line_windows, left_window, right_window):
    left = band_mean_log(x, wl, *left_window)
    right = band_mean_log(x, wl, *right_window)
    line = np.column_stack([band_mean_log(x, wl, *window) for window in line_windows])
    centers = [np.mean(w) for w in line_windows]
    c = (np.mean(centers) - np.mean(left_window)) / (np.mean(right_window) - np.mean(left_window))
    baseline = left * (1.0 - c) + right * c
    return np.nanmean(line, axis=1) - baseline


def transition_time(data):
    gas = smooth_gas(data["gas"])
    total = gas.sum(axis=1)
    co_share = gas[:, 1] / np.maximum(total, 0.1)
    ok = (co_share >= 0.5) & (total > 0.5)
    dwell = int(10 / DT)
    for i in range(len(ok) - dwell + 1):
        if ok[i:i + dwell].all():
            return i * DT
    return np.nan


def wilcoxon_p(values):
    values = np.asarray(values, dtype=float)
    if np.allclose(values, 0):
        return 1.0
    return float(wilcoxon(values).pvalue)


def holm_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * p[idx])
        running = max(running, value)
        adj[idx] = running
    return adj


def spectral_event_diagnostics(heats, output_dir):
    feature_deltas = {name: [] for name in ("Kdoublet", "K766", "K770", "OI777", "band600", "band475", "brightness")}
    controlled_deltas = {name: [] for name in ("Kdoublet", "OI777", "band600")}
    for heat, data in heats.items():
        wl, x = data["wavelength"], data["spectrum"]
        t = np.arange(len(x)) * DT
        tg = transition_time(data)
        if not np.isfinite(tg):
            continue
        f = {
            "K766": local_contrast(x, wl, [(766.0, 767.0)], (762.5, 764.5), (768.0, 768.8)),
            "K770": local_contrast(x, wl, [(769.4, 770.4)], (768.0, 768.8), (772.0, 774.0)),
            "Kdoublet": local_contrast(x, wl, [(765.8, 767.1), (769.3, 770.5)], (762.5, 764.5), (772.0, 774.0)),
            "OI777": local_contrast(x, wl, [(776.2, 778.5)], (772.0, 774.0), (780.0, 782.0)),
            "band600": local_contrast(x, wl, [(600, 605)], (595, 598), (608, 611)),
            "band475": local_contrast(x, wl, [(475, 480)], (470, 473), (482, 485)),
            "brightness": band_mean_log(x, wl, 345, 1015),
        }
        for name, values in f.items():
            values = median_filter(np.nan_to_num(values, nan=np.nanmedian(values)), size=11, mode="nearest")
            pre = values[(t >= tg - 60) & (t < tg - 30)]
            post = values[(t >= tg - 30) & (t < tg)]
            feature_deltas[name].append(float(np.mean(post) - np.mean(pre)))
        # Per-heat trend control: elapsed time, cumulative oxygen, and broad brightness.
        for name in controlled_deltas:
            values = f[name]
            ok = np.isfinite(values)
            design = np.column_stack([np.ones(len(t)), t / 600.0,
                                      (data["oxygen"] - data["oxygen"][0]) / 3000.0,
                                      f["brightness"]])
            coef = np.linalg.lstsq(design[ok], values[ok], rcond=None)[0]
            residual = values - design @ coef
            residual = median_filter(np.nan_to_num(residual, nan=np.nanmedian(residual)), size=11, mode="nearest")
            pre = residual[(t >= tg - 60) & (t < tg - 30)]
            post = residual[(t >= tg - 30) & (t < tg)]
            controlled_deltas[name].append(float(np.mean(post) - np.mean(pre)))

    rows = []
    for name, values in feature_deltas.items():
        rows.append([name, len(values), int(np.sum(np.asarray(values) > 0)), float(np.mean(values)),
                     float(np.median(values)), wilcoxon_p(values), np.nan, np.nan, np.nan, np.nan])
    control_names = list(controlled_deltas)
    control_ps = [wilcoxon_p(controlled_deltas[name]) for name in control_names]
    control_holm = holm_adjust(control_ps)
    for i, name in enumerate(control_names):
        values = controlled_deltas[name]
        idx = next(j for j, row in enumerate(rows) if row[0] == name)
        rows[idx][6] = int(np.sum(np.asarray(values) > 0))
        rows[idx][7] = float(np.mean(values))
        rows[idx][8] = float(control_ps[i])
        rows[idx][9] = float(control_holm[i])
    pd.DataFrame(rows, columns=["feature", "n_heats", "positive_heats", "mean_delta", "median_delta",
                                "raw_wilcoxon_p", "controlled_positive_heats", "controlled_mean_delta",
                                "controlled_wilcoxon_p", "controlled_holm_p"]).to_csv(output_dir / "spectral_event_features.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-zip", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(args.data_zip.read_bytes()).hexdigest()
    heats = load_data_zip(args.data_zip)
    rows = sum(len(v["gas"]) for v in heats.values())
    tail_diagnostics(heats, args.output_dir)
    bounds = stage_diagnostics(heats, args.output_dir)
    stage_lag_diagnostics(heats, bounds, args.output_dir)
    spectral_event_diagnostics(heats, args.output_dir)
    manifest = {
        "data_zip": str(args.data_zip),
        "sha256": digest,
        "heats": len(heats),
        "rows": rows,
        "wavelength_channels": 2047,
        "dt_s": DT,
        "outputs": [
            "tail_dynamics.csv", "tail_window_sensitivity.csv", "stage_boundaries.csv",
            "stage_slopes.csv", "stage_summary.csv", "stage_lags_by_heat.csv",
            "stage_lag_summary.csv", "spectral_event_features.csv",
        ],
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
