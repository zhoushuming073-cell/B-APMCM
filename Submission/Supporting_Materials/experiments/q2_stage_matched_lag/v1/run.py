#!/usr/bin/env python3
"""Q2.3 third-round experiment: within-heat, oxygen-progress-matched lag contrasts.

This is a supplementary, exploratory analysis. It does not overwrite the frozen Q2.2/Q2.3
results. Run from repository root after Q2.1 finalization:
    python experiments/q2_stage_matched_lag/v1/run.py
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest

HEATS = ["6B03205","6B03206","6B03270","6B03271","6C02057",
         "6C02058","6C02059","6C02060","6C02107","6C02108"]
SCENARIOS = ["H0","H1"]
GASES = ["CO2","CO"]
WINDOWS = [90,180]
METRICS = ["pearson","spearman"]
QUANTILE_PAIRS = [(0.25,0.75),(0.30,0.70),(0.40,0.60)]
CALIPERS = [0.05,0.10,0.15,0.20]

def corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x-x.mean(); y=y-y.mean()
    den=np.sqrt((x@x)*(y@y))
    return float(x@y/den) if den>0 else np.nan

def assoc(x,y,metric):
    if metric=="pearson":
        return corr(x,y)
    return corr(pd.Series(x).rank(method="average"),pd.Series(y).rank(method="average"))

def detrend(x):
    x=np.asarray(x,float); t=np.arange(len(x),dtype=float)
    X=np.column_stack([np.ones(len(x)),t])
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

def whiten(x,y):
    x=detrend(x); y=detrend(y)
    den=x[:-1]@x[:-1]
    phi=np.clip((x[:-1]@x[1:])/den if den>0 else 0.0,-0.98,0.98)
    return x[1:]-phi*x[:-1], y[1:]-phi*y[:-1], float(phi)

def best_lag(x,y,maxlag,metric):
    vals=[abs(assoc(x,y[k:k+len(x)],metric)) for k in range(maxlag+1)]
    k=int(np.nanargmax(vals))
    return k,float(vals[k]),float(assoc(x,y[k:k+len(x)],metric))

def estimate_windows(df):
    rows=[]
    for metric in METRICS:
      for scenario in SCENARIOS:
       for window in WINDOWS:
        n_input=window//3
        maxlag=window//6
        for heat in HEATS:
         f=df[(df.scenario==scenario)&(df.heat==heat)].sort_values("time_s")
         if f.empty: raise ValueError(f"Missing heat/scenario: {scenario}/{heat}")
         t=f.time_s.to_numpy(float)
         q=f.oxygen_flow_analytic_m3_s.to_numpy(float)
         O=f.oxygen_cumulative_fit_m3.to_numpy(float)
         for start in np.arange(0,t[-1]-window+1e-8,window):
          grid=start+2*np.arange(window//2)
          qx=np.interp(grid[:n_input],t,q)
          progress=(np.interp(start+window/2,t,O)-O[0])/(O[-1]-O[0])
          dose=float(np.interp(start+window/2,t,O)-O[0])
          for gas in GASES:
           gy=np.interp(grid,t,f[f"{gas}_flow_k2_m3_s"].to_numpy(float))
           if np.std(qx)<1e-7 or np.std(gy)<1e-7: continue
           x,y,phi=whiten(qx,gy)
           if np.std(x)<1e-9 or np.std(y)<1e-9: continue
           k,peak,signed=best_lag(x,y,maxlag,metric)
           rows.append({
             "metric":metric,"scenario":scenario,"window_s":window,"heat":heat,"gas":gas,
             "start_s":float(start),"flow_m3_s":float(qx.mean()),
             "oxygen_progress":float(progress),"oxygen_dose_m3":dose,
             "lag_s":float(2*k),"peak_abs_assoc":peak,"signed_assoc_at_peak":signed,
             "ar1_phi_from_oxygen":phi,"boundary_peak":bool(k==maxlag)
           })
    return pd.DataFrame(rows)

def match_heat(g, q_low, q_high, caliper):
    """Strict non-overlapping low/high tails + one-to-one min-distance matching.

    Ties at quantile cut points are excluded. This prevents a constant-flow window
    from appearing in both treatment and control groups.
    """
    lo=g[g.flow_m3_s<q_low].sort_values("oxygen_progress").reset_index(drop=True)
    hi=g[g.flow_m3_s>q_high].sort_values("oxygen_progress").reset_index(drop=True)
    if lo.empty or hi.empty:
        return []
    dist=np.abs(hi.oxygen_progress.to_numpy()[:,None]-lo.oxygen_progress.to_numpy()[None,:])
    nh,nl=dist.shape
    cost=np.full((nh,nl+nh),caliper+1e-8)
    cost[:,:nl]=np.where(dist<=caliper,dist,caliper+1.0)
    for i in range(nh):
        cost[i,nl+i]=caliper+1e-8
    rr,cc=linear_sum_assignment(cost)
    pairs=[]
    for i,j in zip(rr,cc):
        if j>=nl or dist[i,j]>caliper:
            continue
        h=hi.iloc[i]; l=lo.iloc[j]
        pairs.append({
          "heat":h.heat,"gas":h.gas,"scenario":h.scenario,"window_s":int(h.window_s),
          "metric":h.metric,"high_start_s":h.start_s,"low_start_s":l.start_s,
          "high_progress":h.oxygen_progress,"low_progress":l.oxygen_progress,
          "progress_gap":abs(h.oxygen_progress-l.oxygen_progress),
          "high_flow_m3_s":h.flow_m3_s,"low_flow_m3_s":l.flow_m3_s,
          "flow_diff_m3_s":h.flow_m3_s-l.flow_m3_s,
          "high_lag_s":h.lag_s,"low_lag_s":l.lag_s,
          "lag_diff_high_minus_low_s":h.lag_s-l.lag_s,
          "high_peak_abs_assoc":h.peak_abs_assoc,"low_peak_abs_assoc":l.peak_abs_assoc,
          "high_signed_assoc":h.signed_assoc_at_peak,"low_signed_assoc":l.signed_assoc_at_peak,
          "high_boundary_peak":h.boundary_peak,"low_boundary_peak":l.boundary_peak
        })
    return pairs

def cluster_bootstrap(values, stat, seed, B=20000):
    x=np.asarray(values,float); x=x[np.isfinite(x)]
    if len(x)==0:return (np.nan,np.nan)
    rng=np.random.default_rng(seed)
    idx=rng.integers(0,len(x),size=(B,len(x)))
    z=x[idx]
    sims=np.median(z,axis=1) if stat=="median" else np.mean(z,axis=1)
    return tuple(np.quantile(sims,[0.025,0.975]))

def summarize_matches(pairs, scenario, window, metric, quantiles, caliper):
    rows=[]; heat_rows=[]
    for gas in GASES:
      g=pairs[pairs.gas==gas].copy() if len(pairs) else pd.DataFrame()
      if g.empty:
        rows.append({
          "scenario":scenario,"window_s":window,"metric":metric,
          "q_low":quantiles[0],"q_high":quantiles[1],"caliper_progress":caliper,
          "gas":gas,"n_pairs":0,"n_heats":0
        })
        continue
      for heat,h in g.groupby("heat"):
        heat_rows.append({
          "scenario":scenario,"window_s":window,"metric":metric,
          "q_low":quantiles[0],"q_high":quantiles[1],"caliper_progress":caliper,
          "gas":gas,"heat":heat,"n_pairs":len(h),
          "heat_median_lag_diff_s":float(h.lag_diff_high_minus_low_s.median()),
          "heat_mean_lag_diff_s":float(h.lag_diff_high_minus_low_s.mean()),
          "heat_median_flow_diff_m3_s":float(h.flow_diff_m3_s.median()),
          "heat_median_progress_gap":float(h.progress_gap.median())
        })
      hd=pd.DataFrame([z for z in heat_rows if z["gas"]==gas and
          z["scenario"]==scenario and z["window_s"]==window and
          z["metric"]==metric and z["q_low"]==quantiles[0] and
          z["q_high"]==quantiles[1] and z["caliper_progress"]==caliper])
      effects=hd.heat_median_lag_diff_s.to_numpy(float)
      neg=int(np.sum(effects<0)); nonzero=int(np.sum(effects!=0))
      p_one=binomtest(neg,nonzero,0.5,alternative="greater").pvalue if nonzero else np.nan
      p_two=binomtest(neg,nonzero,0.5,alternative="two-sided").pvalue if nonzero else np.nan
      seed=20261128+int(caliper*1000)+int(quantiles[0]*100)+window+(0 if scenario=="H0" else 100)+(0 if gas=="CO2" else 10)+(0 if metric=="pearson" else 1)
      ci_med=cluster_bootstrap(effects,"median",seed)
      ci_mean=cluster_bootstrap(effects,"mean",seed+1)
      rows.append({
        "scenario":scenario,"window_s":window,"metric":metric,
        "q_low":quantiles[0],"q_high":quantiles[1],"caliper_progress":caliper,
        "gas":gas,"n_pairs":len(g),"n_heats":len(hd),
        "median_progress_gap":float(g.progress_gap.median()),
        "max_progress_gap":float(g.progress_gap.max()),
        "median_flow_diff_m3_s":float(g.flow_diff_m3_s.median()),
        "median_pair_lag_diff_s":float(g.lag_diff_high_minus_low_s.median()),
        "median_heat_lag_diff_s":float(np.median(effects)),
        "mean_heat_lag_diff_s":float(np.mean(effects)),
        "cluster_boot_median_low_s":ci_med[0],"cluster_boot_median_high_s":ci_med[1],
        "cluster_boot_mean_low_s":ci_mean[0],"cluster_boot_mean_high_s":ci_mean[1],
        "heats_high_flow_shorter":neg,
        "heats_high_flow_longer":int(np.sum(effects>0)),
        "heats_equal":int(np.sum(effects==0)),
        "sign_test_p_high_flow_shorter_one_sided":p_one,
        "sign_test_p_two_sided":p_two
      })
    return rows,heat_rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,default=Path("results/q2/final_freeze_2026-11-28/final_curves.csv"))
    ap.add_argument("--output",type=Path,default=Path("results/q2/stage_matched_lag_2026-11-28"))
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    raw=args.input.read_bytes()
    df=pd.read_csv(args.input)
    required={"scenario","heat","time_s","oxygen_cumulative_fit_m3",
              "oxygen_flow_analytic_m3_s","CO2_flow_k2_m3_s","CO_flow_k2_m3_s"}
    missing=required-set(df.columns)
    if missing: raise ValueError(f"Missing columns: {sorted(missing)}")
    if len(df)!=27178 or df.heat.nunique()!=10 or set(df.scenario)!=set(SCENARIOS):
        raise ValueError("Expected 27,178 frozen Q2.1 rows, 10 heats, H0/H1.")
    if df[list(required)].isna().any().any():
        raise ValueError("Missing values in required frozen-curve columns.")

    windows=estimate_windows(df)
    all_pairs=[]; all_summaries=[]; all_heat=[]
    for scenario in SCENARIOS:
      for window in WINDOWS:
       for metric in METRICS:
        sub=windows[(windows.scenario==scenario)&(windows.window_s==window)&(windows.metric==metric)]
        quantile_specs=QUANTILE_PAIRS if window==90 else [(0.25,0.75)]
        for quantiles in quantile_specs:
          for caliper in CALIPERS:
            pairs=[]
            for (heat,gas),g in sub.groupby(["heat","gas"]):
              qlow,qhigh=np.quantile(g.flow_m3_s,[quantiles[0],quantiles[1]])
              for z in match_heat(g,qlow,qhigh,caliper):
                z.update({"q_low":quantiles[0],"q_high":quantiles[1],"caliper_progress":caliper})
                pairs.append(z)
            pair_df=pd.DataFrame(pairs)
            if len(pair_df): all_pairs.append(pair_df)
            s,h=summarize_matches(pair_df,scenario,window,metric,quantiles,caliper)
            all_summaries.extend(s); all_heat.extend(h)

    pair_table=pd.concat(all_pairs,ignore_index=True) if all_pairs else pd.DataFrame()
    summary=pd.DataFrame(all_summaries)
    heat_effects=pd.DataFrame(all_heat)
    windows.to_csv(args.output/"lag_windows.csv",index=False,float_format="%.8g")
    pair_table.to_csv(args.output/"matched_pairs.csv",index=False,float_format="%.8g")
    summary.to_csv(args.output/"matching_summary.csv",index=False,float_format="%.8g")
    heat_effects.to_csv(args.output/"heat_level_effects.csv",index=False,float_format="%.8g")
    manifest={
      "experiment":"Q2.3 third-round within-heat oxygen-progress matching",
      "input_curves_sha256":hashlib.sha256(raw).hexdigest(),
      "input_rows":len(df),"heats":HEATS,"scenarios":SCENARIOS,
      "windows_s":WINDOWS,"metrics":METRICS,"grid_s":2,
      "lag_estimator":"same as frozen Q2 relationship study: first 2/3 window oxygen flow, full-window gas flow, linear detrend, AR(1) prewhitening using oxygen-derived phi, max absolute association over nonnegative lags",
      "matching":"within heat and gas; bottom/top flow quantile tails; strict inequalities exclude tied cutpoints; one-to-one Hungarian minimum-distance matching on cumulative-oxygen progress; no replacement",
      "quantile_pairs":QUANTILE_PAIRS,"calipers_progress":CALIPERS,
      "effect_definition":"high-flow matched-window lag minus low-flow matched-window lag; negative supports shorter lag at high flow",
      "inference":"heat-level median matched difference; percentile cluster bootstrap over heats; exact one-sided sign test for negative differences",
      "source_data_zip_sha256":"cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727",
      "frozen_curve_reference_sha256":"9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49",
      "notes":["Exploratory post-freeze supplement; does not modify frozen Q2.2/Q2.3 outputs.",
               "Caliper 0.15 is a relaxed diagnostic tolerance, not exact stage matching.",
               "Windows are not independent; inference resamples heat clusters.",
               "No matches means lack of common support, not evidence of zero effect."]
    }
    (args.output/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("Wrote:",args.output)
    print(summary[(summary.window_s==90)&(summary.metric=="pearson")&
                  (summary.q_low==0.25)&(summary.q_high==0.75)&
                  (summary.caliper_progress.isin([0.10,0.15]))].to_string(index=False))

if __name__=="__main__":
    main()
