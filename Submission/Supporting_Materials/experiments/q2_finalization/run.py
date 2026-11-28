from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from numpy.lib.stride_tricks import sliding_window_view

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.q2_instantaneous_flow.v4.run import run as run_flow
from experiments.q2_anomaly_review.v1.run import run as run_anomaly


def isolated_spike_flags(t, y):
    """Legacy isolated-spike screen; returns candidates, never edits raw values."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    flags = np.zeros(len(y), dtype=bool)
    if len(y) < 9:
        return flags
    w = sliding_window_view(y, 9)
    fl = w[:, [0, 1, 2, 3, 5, 6, 7, 8]]
    left, right = np.median(fl[:, :4], axis=1), np.median(fl[:, 4:], axis=1)
    med = np.median(fl, axis=1)
    limit = np.maximum(.2, 6 * 1.4826 * np.median(np.abs(fl - med[:, None]), axis=1))
    center = w[:, 4]
    distance = np.minimum.reduce([np.abs(center-w[:,3]), np.abs(center-w[:,5]),
                                  np.abs(center-left), np.abs(center-right)])
    dt = float(np.median(np.diff(t)))
    flags[4:-4] = ((t[8:]-t[:-8] <= 8*dt+1e-8) & (np.abs(left-right)<limit)
                   & (np.abs(w[:,3]-w[:,5])<limit)
                   & ((center-w[:,3])*(center-w[:,5])>0) & (distance>limit))
    return flags


def fit_k2(t, y):
    """2-second local-window means at 2-second knots, followed by PCHIP."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    knots = np.unique(np.r_[t[0], np.arange(t[0]+2.0, t[-1], 2.0), t[-1]])
    centers, means = [t[0]], [y[0]]
    for knot in knots[1:-1]:
        take = np.abs(t-knot) <= 1.0
        if take.sum() >= 2:
            centers.append(float(t[take].mean()))
            means.append(float(y[take].mean()))
    centers.append(t[-1]); means.append(y[-1])
    return PchipInterpolator(centers, np.maximum(means, 0), extrapolate=False)


def run(input_path: Path, out: Path):
    source = pd.read_csv(input_path).sort_values(["attachment", "row"])
    required = {"row", "O", "CO2", "CO", "attachment", "heat"}
    if required-set(source.columns) or source.heat.nunique()!=10 or len(source)!=13589:
        raise ValueError("Expected the corrected 10-heat, 13,589-row process dataset.")

    flow_dir, anomaly_dir = out/"instantaneous_flow", out/"anomaly_review"
    run_flow(input_path, flow_dir)
    run_anomaly(input_path, anomaly_dir)
    flow = pd.read_csv(flow_dir/"flow_profile.csv")
    flow_sum = pd.read_csv(flow_dir/"flow_summary.csv")
    flow_qa = pd.read_csv(flow_dir/"flow_validation_checks.csv")
    review = pd.read_csv(anomaly_dir/"actual_review_queue.csv")

    curves, summaries, sensitivity = [], [], []
    for scenario in ("H0", "H1"):
        for heat, d in source.groupby("heat", sort=False):
            d = d.sort_values("row")
            dt = 1.0 if scenario=="H1" and heat=="6C02057" else .5
            t = (d.row.to_numpy(float)-1)*dt
            oxy = flow[(flow.scenario==scenario)&(flow.heat==heat)].sort_values("row")
            gas_results, gas_stats = {}, {}
            for gas in ("CO2", "CO"):
                raw = d[gas].to_numpy(float)
                primary = fit_k2(t, raw)  # primary curve preserves all raw observations
                screened = ~isolated_spike_flags(t, raw)
                clean = fit_k2(t[screened], raw[screened])
                yhat = primary(t)
                anti = primary.antiderivative()
                cumulative = anti(t)-anti(t[0])
                clean_anti = clean.antiderivative()
                raw_total = float(trapezoid(raw, x=t))
                total = float(cumulative[-1])
                gas_results[gas] = (yhat, cumulative)
                gas_stats[gas] = {
                    "raw_total": raw_total, "k2_total": total,
                    "difference_pct": (total/raw_total-1)*100,
                    "flagged_points": int((~screened).sum())
                }
                sensitivity.append({
                    "scenario": scenario, "heat": heat, "gas": gas,
                    "isolated_candidate_points": int((~screened).sum()),
                    "max_abs_curve_difference_m3_s": float(np.max(np.abs(yhat-clean(t)))),
                    "cumulative_integral_difference_m3": float(total-(clean_anti(t[-1])-clean_anti(t[0]))),
                    "primary_uses_raw_points": True
                })

            curves.append(pd.DataFrame({
                "scenario":scenario, "heat":heat, "row":d.row.to_numpy(int), "time_s":t,
                "oxygen_cumulative_raw_m3":d.O.to_numpy(float),
                "oxygen_cumulative_fit_m3":oxy.cumulative_O_fit_m3.to_numpy(float),
                "oxygen_flow_analytic_m3_s":oxy.flow_analytic_best_m3_s.to_numpy(float),
                "CO2_flow_raw_m3_s":d.CO2.to_numpy(float),
                "CO2_flow_k2_m3_s":gas_results["CO2"][0],
                "CO2_cumulative_k2_m3":gas_results["CO2"][1],
                "CO_flow_raw_m3_s":d.CO.to_numpy(float),
                "CO_flow_k2_m3_s":gas_results["CO"][0],
                "CO_cumulative_k2_m3":gas_results["CO"][1]
            }))

            fs = flow_sum[(flow_sum.scenario==scenario)&(flow_sum.heat==heat)].iloc[0]
            fq = flow_qa[(flow_qa.scenario==scenario)&(flow_qa.heat==heat)].iloc[0]
            summaries.append({
                "scenario":scenario, "heat":heat, "dt_s_per_row":dt, "n_rows":len(d), "duration_s":t[-1],
                "oxygen_knot_spacing_s":fs.selected_spacing_s, "oxygen_penalty":fs.selected_penalty,
                "oxygen_cumulative_CV_MAE_m3":fs.full_CV_MAE_m3,
                "oxygen_flow_min_m3_s":fs.flow_min_m3_s, "oxygen_flow_median_m3_s":fs.flow_median_m3_s,
                "oxygen_flow_mean_m3_s":fs.flow_mean_m3_s, "oxygen_flow_p95_m3_s":fs.flow_p95_m3_s,
                "oxygen_flow_max_m3_s":fs.flow_max_m3_s,
                "oxygen_mass_balance_error_m3":fq.abs_numeric_mass_balance_error_m3,
                "oxygen_nonnegative_pass":fq.pass_nonnegative, "oxygen_monotonic_pass":fq.pass_monotonic,
                "CO2_raw_total_m3":gas_stats["CO2"]["raw_total"], "CO2_k2_total_m3":gas_stats["CO2"]["k2_total"],
                "CO2_k2_difference_pct":gas_stats["CO2"]["difference_pct"],
                "CO2_isolated_candidate_points":gas_stats["CO2"]["flagged_points"],
                "CO_raw_total_m3":gas_stats["CO"]["raw_total"], "CO_k2_total_m3":gas_stats["CO"]["k2_total"],
                "CO_k2_difference_pct":gas_stats["CO"]["difference_pct"],
                "CO_isolated_candidate_points":gas_stats["CO"]["flagged_points"],
                "anomaly_review_events_for_heat":int((review.heat==heat).sum())
            })

    final_curves = pd.concat(curves, ignore_index=True)
    final_summary = pd.DataFrame(summaries)
    final_curves.to_csv(out/"final_curves.csv", index=False)
    final_summary.to_csv(out/"final_heat_summary.csv", index=False)
    flow_qa.to_csv(out/"flow_validation_checks.csv", index=False)
    pd.DataFrame(sensitivity).to_csv(out/"gas_k2_raw_vs_flagged_check.csv", index=False)

    counts = review.groupby("anomaly_type").size().to_dict()
    disposition = [{"anomaly_type":k,"event_count":int(v),
                    "decision":"pending_manual_review; preserve_raw; no_automatic_correction"}
                   for k,v in sorted(counts.items())]
    disposition.append({"anomaly_type":"TOTAL","event_count":int(len(review)),
                        "decision":"all_candidates_pending; raw_data_unchanged"})
    pd.DataFrame(disposition).to_csv(out/"anomaly_disposition_summary.csv", index=False)

    manifest = {
        "process_csv_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "heats":10,"input_rows":len(source),"timebase_scenarios":["H0","H1"],
        "final_curve_rows":len(final_curves),"anomaly_candidate_events":len(review),
        "raw_values_modified":False,"automatic_anomaly_corrections":0
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")

    assert len(final_curves)==27178 and len(final_summary)==20
    assert flow_qa.fit_success.all() and flow_qa.pass_nonnegative.all()
    assert flow_qa.pass_monotonic.all() and flow_qa.pass_mass_balance_0_1m3.all()
    assert not final_curves.isna().any().any()
    print("Q2.1 FINAL FREEZE PASSED")
    print(f"Input rows: {len(source):,}; final curve rows: {len(final_curves):,}")
    print(f"Flow QA: {len(flow_qa)}/20 passed; max mass-balance error: {flow_qa.abs_numeric_mass_balance_error_m3.max():.6f} m3")
    print(f"Anomaly candidates: {len(review)}; raw data changed: no")
    print(f"Output: {out}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Reproduce frozen Q2.1 six-curve outputs and QA")
    parser.add_argument("--input",type=Path,default=ROOT/"data/interim/extracted/process.csv")
    parser.add_argument("--out",type=Path,default=ROOT/"results/q2/final_freeze_2026-11-28")
    args=parser.parse_args()
    run(args.input,args.out)
