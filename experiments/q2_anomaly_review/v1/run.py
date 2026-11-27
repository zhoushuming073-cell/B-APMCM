from __future__ import annotations

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter


def isolated_spike_flags(t, values):
    """Legacy isolated one-point spike detector, retained unchanged in logic."""
    v = np.asarray(values, dtype=float)
    flag = np.zeros(len(v), dtype=bool)
    if len(v) < 9:
        return flag
    windows = np.lib.stride_tricks.sliding_window_view(v, 9)
    flanks = windows[:, [0, 1, 2, 3, 5, 6, 7, 8]]
    left = np.median(flanks[:, :4], axis=1)
    right = np.median(flanks[:, 4:], axis=1)
    med = np.median(flanks, axis=1)
    limit = np.maximum(
        0.2,
        6 * 1.4826 * np.median(np.abs(flanks - med[:, None]), axis=1),
    )
    center = windows[:, 4]
    distance = np.minimum.reduce([
        np.abs(center - windows[:, 3]),
        np.abs(center - windows[:, 5]),
        np.abs(center - left),
        np.abs(center - right),
    ])
    flag[4:-4] = (
        ((t[8:] - t[:-8]) <= 8 * np.median(np.diff(t)) + 1e-8)
        & (np.abs(left - right) < limit)
        & (np.abs(windows[:, 3] - windows[:, 5]) < limit)
        & ((center - windows[:, 3]) * (center - windows[:, 5]) > 0)
        & (distance > limit)
    )
    return flag


def shift_scores(values, window=20):
    """Two-sided local linear extrapolation; returns signed level jump and slope mismatch."""
    z = np.log1p(np.maximum(np.asarray(values, dtype=float), 0))
    n = len(z)
    jump = np.full(n, np.nan)
    slope_gap = np.full(n, np.nan)
    if n <= 2 * window:
        return jump, slope_gap

    cs = np.r_[0.0, np.cumsum(z)]
    idx = np.arange(n, dtype=float)
    csi = np.r_[0.0, np.cumsum(idx * z)]
    i = np.arange(window, n - window)

    pre_sum = cs[i] - cs[i - window]
    pre_ix = csi[i] - csi[i - window]
    post_sum = cs[i + window] - cs[i]
    post_ix = csi[i + window] - csi[i]

    center_pre = i - window + (window - 1) / 2
    center_post = i + (window - 1) / 2
    sxx = window * (window * window - 1) / 12

    slope_pre = (pre_ix - center_pre * pre_sum) / sxx
    slope_post = (post_ix - center_post * post_sum) / sxx
    mean_pre = pre_sum / window
    mean_post = post_sum / window

    pred_pre_at_split = mean_pre + slope_pre * (i - center_pre)
    pred_post_at_split = mean_post + slope_post * (i - center_post)

    jump[i] = pred_post_at_split - pred_pre_at_split
    slope_gap[i] = np.abs(slope_post - slope_pre) * window
    return jump, slope_gap


def contiguous_runs(mask):
    indices = np.flatnonzero(mask)
    if len(indices) == 0:
        return []
    groups = np.split(indices, np.where(np.diff(indices) > 1)[0] + 1)
    return [(int(g[0]), int(g[-1])) for g in groups if len(g)]


def detect_series(values, dt_s, sigma_override=None):
    v = np.asarray(values, dtype=float)
    n = len(v)
    z = np.log1p(np.maximum(v, 0))
    baseline_log = median_filter(z, size=21, mode="nearest")
    residual = z - baseline_log

    # Quantile scale is used because staircase/quantized measurements often have MAD=0.
    sigma = (max(float(np.quantile(np.abs(residual), 0.90) / 1.645), 0.003)
             if sigma_override is None else max(float(sigma_override), 0.003))
    burst_threshold = max(3.5 * sigma, 0.018)
    shift_threshold = max(8.0 * sigma, 0.040)
    slope_gap_threshold = max(8.0 * sigma, 0.040)

    t = np.arange(n, dtype=float) * dt_s
    flags = {
        "isolated_spike": isolated_spike_flags(t, v),
        "multi_point_burst": np.zeros(n, dtype=bool),
        "sustained_shift": np.zeros(n, dtype=bool),
    }

    burst_events = []
    for start, end in contiguous_runs(np.abs(residual) > burst_threshold):
        length = end - start + 1
        if not 2 <= length <= 8:
            continue
        segment = residual[start:end + 1]
        sign = np.sign(np.median(segment))
        if sign == 0 or not np.all(np.sign(segment) == sign):
            continue
        before = residual[max(0, start - 10):start]
        after = residual[end + 1:min(n, end + 11)]
        if len(before) < 6 or len(after) < 6:
            continue
        if abs(np.median(before)) >= burst_threshold * 0.5:
            continue
        if abs(np.median(after)) >= burst_threshold * 0.5:
            continue
        flags["multi_point_burst"][start:end + 1] = True
        burst_events.append((start, end, float(np.median(segment))))

    jump, slope_gap = shift_scores(v, window=20)
    candidates = np.flatnonzero(
        np.isfinite(jump)
        & (np.abs(jump) > shift_threshold)
        & (slope_gap < slope_gap_threshold)
    )
    shift_events = []
    if len(candidates):
        groups = np.split(candidates, np.where(np.diff(candidates) > 10)[0] + 1)
        for group in groups:
            split = int(group[np.argmax(np.abs(jump[group]))])
            flags["sustained_shift"][split] = True
            shift_events.append((
                split, float(jump[split]), float(slope_gap[split])
            ))

    metadata = {
        "sigma_log": sigma,
        "burst_threshold_log": burst_threshold,
        "shift_threshold_log": shift_threshold,
        "slope_gap_threshold_log": slope_gap_threshold,
        "baseline_log": baseline_log,
        "residual_log": residual,
        "burst_events": burst_events,
        "shift_events": shift_events,
    }
    return flags, metadata


def inject_series(values, rng, anomaly_type, amplitude, duration=1, n_events=5):
    x = np.asarray(values, dtype=float).copy()
    occupied = np.zeros(len(x), dtype=bool)
    events = []
    n = len(x)
    for _ in range(n_events):
        for _attempt in range(500):
            lower = max(40, int(0.55 * n))
            if anomaly_type == "sustained_shift":
                start = int(rng.integers(lower, n - 40))
                end = n - 1
            else:
                start = int(rng.integers(lower, n - duration - 40))
                end = start + duration - 1
            if occupied[max(0, start - 25):min(n, end + 26)].any():
                continue
            local_level = float(np.median(x[max(0, start - 5):min(n, start + 5)]))
            if local_level < amplitude + 2:
                continue
            sign = int(rng.choice([-1, 1]))
            if anomaly_type == "sustained_shift":
                x[start:] += sign * amplitude
            else:
                x[start:end + 1] += sign * amplitude
            x = np.maximum(x, 0)
            occupied[start:end + 1] = True
            events.append((start, end, sign, amplitude))
            break
    return x, events


def build_actual_outputs(frame, output_dir):
    review_rows, point_rows, threshold_rows = [], [], []
    for heat, d in frame.groupby("heat", sort=False):
        d = d.sort_values("row")
        n = len(d)
        dt_h0 = 0.5
        dt_h1 = 1.0 if heat == "6C02057" else 0.5
        for gas in ("CO2", "CO"):
            values = d[gas].to_numpy(float)
            flags, meta = detect_series(values, dt_h0)
            threshold_rows.append({
                "heat": heat,
                "gas": gas,
                "sigma_log": meta["sigma_log"],
                "burst_threshold_log": meta["burst_threshold_log"],
                "shift_threshold_log": meta["shift_threshold_log"],
                "slope_gap_threshold_log": meta["slope_gap_threshold_log"],
                "isolated_spike_points": int(flags["isolated_spike"].sum()),
                "burst_points": int(flags["multi_point_burst"].sum()),
                "shift_boundaries": int(flags["sustained_shift"].sum()),
                "timebase_note": "row-space detector; H0/H1 changes displayed time only",
            })

            # Event-level queue. Every candidate remains pending manual review.
            for i in np.flatnonzero(flags["isolated_spike"]):
                before = float(np.median(values[max(0, i - 4):i])) if i >= 1 else np.nan
                after = float(np.median(values[i + 1:min(n, i + 5)])) if i + 1 < n else np.nan
                review_rows.append({
                    "heat": heat, "gas": gas, "anomaly_type": "isolated_spike",
                    "event_start_row": int(d["row"].iloc[i]),
                    "event_end_row": int(d["row"].iloc[i]),
                    "review_window_end_row": int(d["row"].iloc[i]),
                    "time_start_H0_s": float(i * dt_h0),
                    "time_start_H1_s": float(i * dt_h1),
                    "observed_median_m3_s": float(values[i]),
                    "baseline_reference_m3_s": float(np.nanmean([before, after])),
                    "signed_log_deviation": float(meta["residual_log"][i]),
                    "score_ratio": np.nan,
                    "review_priority": "medium",
                    "review_status": "pending",
                    "automatic_action": "legacy isolated-spike filter may exclude from cleaned reconstruction; raw value retained",
                    "review_note": "Check neighbouring readings and process log before treating as sensor error.",
                })

            for start, end, amplitude_log in meta["burst_events"]:
                baseline = float(np.expm1(np.median(meta["baseline_log"][start:end + 1])))
                review_rows.append({
                    "heat": heat, "gas": gas, "anomaly_type": "multi_point_burst",
                    "event_start_row": int(d["row"].iloc[start]),
                    "event_end_row": int(d["row"].iloc[end]),
                    "review_window_end_row": int(d["row"].iloc[end]),
                    "time_start_H0_s": float(start * dt_h0),
                    "time_start_H1_s": float(start * dt_h1),
                    "observed_median_m3_s": float(np.median(values[start:end + 1])),
                    "baseline_reference_m3_s": baseline,
                    "signed_log_deviation": amplitude_log,
                    "score_ratio": abs(amplitude_log) / meta["burst_threshold_log"],
                    "review_priority": "high" if abs(amplitude_log) >= 0.05 else "medium",
                    "review_status": "pending",
                    "automatic_action": "flag only; no automatic correction",
                    "review_note": "Could be a short sensor disturbance or a real process transient; retain until reviewed.",
                })

            for split, jump, slope_gap in meta["shift_events"]:
                before = float(np.median(values[max(0, split - 20):split]))
                after = float(np.median(values[split:min(n, split + 20)]))
                late_collapse = (
                    split >= int(0.90 * n) and jump < -0.50 and after < 1.0
                )
                shift_note = (
                    "Late-series collapse to near zero: possible process stop/end-stage transition or sensor offset. "
                    "Check operating log before any correction."
                    if late_collapse else
                    "Possible sensor offset or genuine process-stage transition; inspect operating log before correction."
                )
                review_rows.append({
                    "heat": heat, "gas": gas, "anomaly_type": "sustained_level_shift",
                    "event_start_row": int(d["row"].iloc[split]),
                    "event_end_row": int(d["row"].iloc[split]),
                    "review_window_end_row": int(d["row"].iloc[min(n - 1, split + 19)]),
                    "time_start_H0_s": float(split * dt_h0),
                    "time_start_H1_s": float(split * dt_h1),
                    "observed_median_m3_s": after,
                    "baseline_reference_m3_s": before,
                    "signed_log_deviation": jump,
                    "score_ratio": abs(jump) / meta["shift_threshold_log"],
                    "slope_gap_log_per_window": slope_gap,
                    "review_priority": "high" if abs(jump) >= 0.15 else "medium",
                    "review_status": "pending",
                    "automatic_action": "flag only; no automatic correction",
                    "review_note": shift_note,
                })

            combined = np.zeros(n, dtype=bool)
            for arr in flags.values():
                combined |= arr
            for i in np.flatnonzero(combined):
                kinds = [name for name, arr in flags.items() if arr[i]]
                point_rows.append({
                    "heat": heat,
                    "gas": gas,
                    "row": int(d["row"].iloc[i]),
                    "time_H0_s": float(i * dt_h0),
                    "time_H1_s": float(i * dt_h1),
                    "observed_median_m3_s": float(values[i]),
                    "local_median_reference_m3_s": float(np.expm1(meta["baseline_log"][i])),
                    "residual_log": float(meta["residual_log"][i]),
                    "flag_types": ";".join(kinds),
                    "manual_review_required": True,
                    "raw_value_preserved": True,
                    "auto_corrected": False,
                })

    review = pd.DataFrame(review_rows).sort_values(
        ["heat", "gas", "event_start_row", "anomaly_type"]
    )
    points = pd.DataFrame(point_rows).sort_values(["heat", "gas", "row"])
    thresholds = pd.DataFrame(threshold_rows)
    review.to_csv(output_dir / "actual_review_queue.csv", index=False)
    points.to_csv(output_dir / "actual_row_flags.csv", index=False)
    thresholds.to_csv(output_dir / "detector_thresholds.csv", index=False)
    return review, points, thresholds


def run_injection_tests(frame, output_dir):
    rng = np.random.default_rng(20261128)
    rows = []
    for heat, d in frame.groupby("heat", sort=False):
        d = d.sort_values("row")
        dt = 1.0 if heat == "6C02057" else 0.5
        for gas in ("CO2", "CO"):
            values = d[gas].to_numpy(float)
            baseline_flags, baseline_meta = detect_series(values, dt)
            cases = (
                ("isolated_spike", 1, 1.0, 5),
                ("isolated_spike", 1, 4.0, 5),
                ("multi_point_burst", 3, 1.0, 5),
                ("multi_point_burst", 3, 4.0, 5),
                ("multi_point_burst", 5, 1.0, 5),
                ("multi_point_burst", 5, 4.0, 5),
                ("sustained_shift", 40, 1.0, 1),
                ("sustained_shift", 40, 4.0, 1),
            )
            for kind, duration, amplitude, n_events in cases:
                test, events = inject_series(
                    values, rng, kind, amplitude, duration=duration, n_events=n_events
                )
                flags, _ = detect_series(
                    test, dt, sigma_override=baseline_meta["sigma_log"]
                )
                hits = []
                allowed = np.zeros(len(values), dtype=bool)
                for start, end, sign, amp in events:
                    if kind == "isolated_spike":
                        hits.append(bool(flags["isolated_spike"][start]))
                        allowed[max(0, start - 4):min(len(values), end + 5)] = True
                    elif kind == "multi_point_burst":
                        hits.append(bool(flags["multi_point_burst"][start:end + 1].any()))
                        allowed[max(0, start - 10):min(len(values), end + 11)] = True
                    else:
                        boundaries = np.flatnonzero(flags["sustained_shift"])
                        hits.append(bool(np.any(np.abs(boundaries - start) <= 10)))
                        allowed[max(0, start - 10):min(len(values), start + 11)] = True

                extra_flags = 0
                for detector_name, detector_flags in flags.items():
                    extra_flags += int(np.sum(
                        detector_flags & ~allowed & ~baseline_flags[detector_name]
                    ))
                rows.append({
                    "heat": heat,
                    "gas": gas,
                    "anomaly_type": kind,
                    "duration_rows": duration,
                    "amplitude_units": amplitude,
                    "injected_events": len(events),
                    "event_recall": float(np.mean(hits)) if hits else np.nan,
                    "new_flags_outside_event_buffer": extra_flags,
                })

    by_heat = pd.DataFrame(rows)
    summary = by_heat.groupby(
        ["anomaly_type", "duration_rows", "amplitude_units"], as_index=False
    ).agg(
        injected_events=("injected_events", "sum"),
        mean_event_recall=("event_recall", "mean"),
        min_event_recall=("event_recall", "min"),
        max_event_recall=("event_recall", "max"),
        mean_new_flags_outside_event_buffer=("new_flags_outside_event_buffer", "mean"),
        max_new_flags_outside_event_buffer=("new_flags_outside_event_buffer", "max"),
    )
    by_heat.to_csv(output_dir / "injection_test_by_heat.csv", index=False)
    summary.to_csv(output_dir / "injection_test_summary.csv", index=False)
    return by_heat, summary


def run(input_path: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path).sort_values(["attachment", "row"])
    required = {"row", "O", "CO2", "CO", "attachment", "heat"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if frame["heat"].nunique() != 10:
        raise ValueError("Expected ten heats in the process dataset.")

    review, points, thresholds = build_actual_outputs(frame, output_dir)
    _, injection_summary = run_injection_tests(frame, output_dir)

    print("Actual review candidates:")
    print(review.groupby(["anomaly_type", "gas"]).size().to_string())
    print("\nSynthetic injection tests:")
    print(injection_summary.round(3).to_string(index=False))
    print("\nNo multi-point or sustained-shift candidate is automatically deleted.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[3]
    parser.add_argument("--input", type=Path, default=root / "data/interim/extracted/process.csv")
    parser.add_argument("--out", type=Path, default=root / "results/q2/anomaly_review/v1")
    args = parser.parse_args()
    run(args.input, args.out)
