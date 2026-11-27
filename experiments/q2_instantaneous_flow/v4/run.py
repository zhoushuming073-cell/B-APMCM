from __future__ import annotations

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import lsq_linear
from scipy.integrate import trapezoid


def fit_integral_bspline(O, t, spacing_s, penalty, train_mask, irls_iter=2):
    """Fit q(t)=B(t)c, c>=0; cumulative O_hat is its analytic integral."""
    O = np.asarray(O, dtype=float)
    t = np.asarray(t, dtype=float)
    train_mask = np.asarray(train_mask, dtype=bool)

    knots = np.r_[
        [t[0]] * 4,
        np.arange(t[0] + spacing_s, t[-1], spacing_s),
        [t[-1]] * 4,
    ]
    n_basis = len(knots) - 4
    basis = BSpline(knots, np.eye(n_basis), 3, extrapolate=False)
    B = np.asarray(basis(t))
    anti = basis.antiderivative()
    A = np.asarray(anti(t) - anti(t[0]))
    D2 = np.diff(np.eye(n_basis), n=2, axis=0)

    weights = np.ones(train_mask.sum())
    coef = np.zeros(n_basis)
    success = False
    for _ in range(irls_iter):
        X = np.vstack([
            A[train_mask] * np.sqrt(weights)[:, None],
            np.sqrt(penalty) * D2,
        ])
        y = np.r_[
            (O[train_mask] - O[0]) * np.sqrt(weights),
            np.zeros(len(D2)),
        ]
        fit = lsq_linear(X, y, bounds=(0, np.inf), max_iter=120)
        coef = fit.x
        success = bool(fit.success)

        residual = O[train_mask] - O[0] - A[train_mask] @ coef
        center = np.median(residual)
        scale = max(1.4826 * np.median(np.abs(residual - center)), 0.1)
        weights = np.minimum(
            1.0,
            1.345 * scale / np.maximum(np.abs(residual - center), 1e-9),
        )

    q = B @ coef                         # analytic derivative dO_hat/dt
    O_hat = O[0] + A @ coef              # exact antiderivative evaluated at t
    return q, O_hat, success


def make_blocks(n):
    blocks = []
    for frac in (0.12, 0.27, 0.42, 0.57, 0.72, 0.87):
        start = max(20, int(n * frac))
        end = min(n - 20, start + 16)
        blocks.append(np.arange(start, end))
    return blocks


def cv_grid(O, t, base_mask, folds, grid):
    rows = []
    for spacing_s, penalty in grid:
        fold_mae = []
        for held in folds:
            train_mask = base_mask.copy()
            train_mask[held] = False
            _, O_hat, _ = fit_integral_bspline(
                O, t, spacing_s, penalty, train_mask, irls_iter=2
            )
            fold_mae.append(float(np.mean(np.abs(O_hat[held] - O[held]))))
        rows.append({
            "spacing_s": spacing_s,
            "penalty": penalty,
            "cv_MAE_m3": float(np.mean(fold_mae)),
            "cv_MAE_se_m3": float(np.std(fold_mae, ddof=1) / np.sqrt(len(fold_mae))),
            **{f"fold_{i+1}_MAE_m3": value for i, value in enumerate(fold_mae)},
        })
    return pd.DataFrame(rows).sort_values(
        ["cv_MAE_m3", "spacing_s", "penalty"]
    ).reset_index(drop=True)


def score(y, pred):
    err = np.asarray(pred) - np.asarray(y)
    denom = np.sum((np.asarray(y) - np.mean(y)) ** 2)
    return {
        "MAE_m3": float(np.mean(np.abs(err))),
        "RMSE_m3": float(np.sqrt(np.mean(err ** 2))),
        "R2": float(1 - np.sum(err ** 2) / denom) if denom > 0 else np.nan,
    }


def run(input_path: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(input_path).sort_values(["attachment", "row"])
    required = {"row", "O", "CO2", "CO", "attachment", "heat"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    heats = list(frame["heat"].drop_duplicates())
    if len(heats) != 10 or heats[4] != "6C02057":
        raise ValueError(f"Unexpected heat order/count: {heats}")

    grid = [(s, lam) for s in (5.0, 10.0, 20.0) for lam in (0.1, 1.0, 10.0, 100.0)]
    profile_rows, summary_rows, sensitivity_rows, cv_rows = [], [], [], []

    for scenario in ("H0", "H1"):
        for heat in heats:
            d = frame[frame["heat"] == heat].sort_values("row")
            O = d["O"].to_numpy(float)
            row = d["row"].to_numpy(float)
            dt = 1.0 if scenario == "H1" and heat == "6C02057" else 0.5
            t = (row - 1.0) * dt
            n = len(O)
            blocks = make_blocks(n)

            # Full-series CV chooses the final reconstruction and a one-SE sensitivity set.
            full_mask = np.ones(n, dtype=bool)
            full_cv_folds = [blocks[k] for k in (0, 2, 3, 5)]
            full_cv = cv_grid(O, t, full_mask, full_cv_folds, grid)
            best = full_cv.iloc[0]
            one_se_cutoff = float(best.cv_MAE_m3 + best.cv_MAE_se_m3)
            near = full_cv[full_cv["cv_MAE_m3"] <= one_se_cutoff].copy()
            if near.empty:
                near = full_cv.iloc[[0]].copy()

            spacing_best = float(best.spacing_s)
            penalty_best = float(best.penalty)
            q_best, O_best, success = fit_integral_bspline(
                O, t, spacing_best, penalty_best, full_mask, irls_iter=2
            )

            q_models = []
            for _, candidate in near.iterrows():
                q_candidate, _, _ = fit_integral_bspline(
                    O, t, float(candidate.spacing_s), float(candidate.penalty),
                    full_mask, irls_iter=2
                )
                q_models.append(q_candidate)
                sensitivity_rows.append({
                    "scenario": scenario,
                    "heat": heat,
                    "spacing_s": float(candidate.spacing_s),
                    "penalty": float(candidate.penalty),
                    "cv_MAE_m3": float(candidate.cv_MAE_m3),
                    "cv_MAE_se_m3": float(candidate.cv_MAE_se_m3),
                    "within_one_SE": True,
                    "is_best": (
                        float(candidate.spacing_s) == spacing_best
                        and float(candidate.penalty) == penalty_best
                    ),
                })
            q_models = np.asarray(q_models)
            q_median = np.median(q_models, axis=0)
            q_p10 = np.quantile(q_models, 0.10, axis=0)
            q_p90 = np.quantile(q_models, 0.90, axis=0)

            def secant_rate(seconds):
                lag = max(1, int(round(seconds / dt)))
                left = np.maximum(0, np.arange(n) - lag)
                return (O - O[left]) / np.maximum(t - t[left], dt)

            q10, q30, q60 = secant_rate(10), secant_rate(30), secant_rate(60)
            for i in range(n):
                profile_rows.append({
                    "scenario": scenario,
                    "heat": heat,
                    "row": int(row[i]),
                    "time_s": float(t[i]),
                    "cumulative_O_raw_m3": float(O[i]),
                    "cumulative_O_fit_m3": float(O_best[i]),
                    "flow_analytic_best_m3_s": float(q_best[i]),
                    "flow_sensitivity_median_m3_s": float(q_median[i]),
                    "flow_sensitivity_p10_m3_s": float(q_p10[i]),
                    "flow_sensitivity_p90_m3_s": float(q_p90[i]),
                    "observed_secant_10s_m3_s": float(q10[i]),
                    "observed_secant_30s_m3_s": float(q30[i]),
                    "observed_secant_60s_m3_s": float(q60[i]),
                })

            integral_q = float(trapezoid(q_best, t))
            summary_rows.append({
                "scenario": scenario,
                "heat": heat,
                "dt_s_per_row": dt,
                "n_rows": n,
                "selected_spacing_s": spacing_best,
                "selected_penalty": penalty_best,
                "full_CV_MAE_m3": float(best.cv_MAE_m3),
                "full_CV_MAE_SE_m3": float(best.cv_MAE_se_m3),
                "near_optimal_model_count": int(len(near)),
                "flow_min_m3_s": float(np.min(q_best)),
                "flow_median_m3_s": float(np.median(q_best)),
                "flow_mean_m3_s": float(np.mean(q_best)),
                "flow_p95_m3_s": float(np.quantile(q_best, 0.95)),
                "flow_max_m3_s": float(np.max(q_best)),
                "integral_q_trapezoid_m3": integral_q,
                "fitted_cumulative_delta_m3": float(O_best[-1] - O_best[0]),
                "numeric_mass_balance_error_m3": integral_q - float(O_best[-1] - O_best[0]),
                "mean_one_SE_sensitivity_band_width_m3_s": float(np.mean(q_p90 - q_p10)),
                "fit_success": success,
            })
            for _, r in full_cv.iterrows():
                cv_rows.append({
                    "scenario": scenario,
                    "heat": heat,
                    "spacing_s": float(r.spacing_s),
                    "penalty": float(r.penalty),
                    "cv_MAE_m3": float(r.cv_MAE_m3),
                    "cv_MAE_se_m3": float(r.cv_MAE_se_m3),
                })

    pd.DataFrame(profile_rows).to_csv(output_dir / "flow_profile.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(output_dir / "flow_summary.csv", index=False)
    pd.DataFrame(sensitivity_rows).to_csv(output_dir / "flow_one_se_sensitivity_set.csv", index=False)
    pd.DataFrame(cv_rows).to_csv(output_dir / "flow_cv_grid.csv", index=False)
    print("\nInstantaneous flow summary:")
    print(pd.DataFrame(summary_rows)[[
        "scenario", "heat", "flow_median_m3_s", "flow_mean_m3_s",
        "flow_p95_m3_s", "mean_one_SE_sensitivity_band_width_m3_s"
    ]].round(4).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[3]
    parser.add_argument("--input", type=Path, default=root / "data/interim/extracted/process.csv")
    parser.add_argument("--out", type=Path, default=root / "results/q2/instantaneous_flow/v4")
    args = parser.parse_args()
    run(args.input, args.out)
