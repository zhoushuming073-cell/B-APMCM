#!/usr/bin/env python3
"""Q2.3 round-2 unified lag-range and estimator-sensitivity experiment."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

HEATS = ["6B03205","6B03206","6B03270","6B03271","6C02057",
         "6C02058","6C02059","6C02060","6C02107","6C02108"]
SCENARIOS = ("H0","H1")
WINDOWS = (90,120,180)
GRID = 2.0
MAX_LAG = 30.0
METRICS = ("pearson","spearman")
PEAKS = ("positive","absolute")
WHITEN = (False,True)

def detrend(x):
    x=np.asarray(x,float); t=np.arange(len(x),dtype=float)
    X=np.column_stack([np.ones(len(x)),t])
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

def corr(x,y,metric):
    x=np.asarray(x,float); y=np.asarray(y,float)
    if metric=="spearman":
        x=rankdata(x,method="average"); y=rankdata(y,method="average")
    x=x-x.mean(); y=y-y.mean()
    den=np.sqrt(float(x@x)*float(y@y))
    return float(x@y/den) if den>0 else np.nan

def prepare(x,y,do_whiten):
    x=detrend(x); y=detrend(y); phi=np.nan
    if do_whiten:
        den=float(x[:-1]@x[:-1])
        phi=float(np.clip((x[:-1]@x[1:])/den if den>0 else 0.,-.98,.98))
        x=x[1:]-phi*x[:-1]
        y=y[1:]-phi*y[:-1]
    return x,y,phi

def cdiff(t,v,half):
    left=np.maximum(t[0],t-half); right=np.minimum(t[-1],t+half)
    return (np.interp(right,t,v)-np.interp(left,t,v))/np.maximum(right-left,1e-12)

def flow_methods(f):
    t=f.time_s.to_numpy(float)
    return {
      "analytic":f.oxygen_flow_analytic_m3_s.to_numpy(float),
      "raw_fd_2s":cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),1.),
      "fit_fd_2s":cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),1.),
      "raw_fd_30s":cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),15.),
      "fit_fd_30s":cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),15.),
      "raw_fd_60s":cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),30.),
      "fit_fd_60s":cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),30.),
    }

def peak(x,y,metric,mode):
    n=len(x); steps=int(MAX_LAG/GRID)
    r=np.array([corr(x,y[k:k+n],metric) for k in range(steps+1)])
    ok=np.isfinite(r)
    if not ok.any(): return np.nan,np.nan,np.nan,False,True
    if mode=="positive":
        ix=np.flatnonzero(ok & (r>0))
        if not len(ix): return np.nan,np.nan,float(np.nanmax(r)),False,True
        k=int(ix[np.argmax(r[ix])])
    else:
        k=int(np.nanargmax(np.abs(r)))
    return k*GRID,float(r[k]),abs(float(r[k])),k==steps,False

def build(curves):
    rows=[]
    for (scen,heat),f in curves.groupby(["scenario","heat"],sort=False):
        f=f.sort_values("time_s"); t=f.time_s.to_numpy(float)
        O=f.oxygen_cumulative_fit_m3.to_numpy(float); span=O[-1]-O[0]
        flows=flow_methods(f)
        gases={"CO":f.CO_flow_k2_m3_s.to_numpy(float),
               "CO2":f.CO2_flow_k2_m3_s.to_numpy(float)}
        for win in WINDOWS:
            grid=np.arange(0.,win,GRID); n=win//3; steps=int(MAX_LAG/GRID)
            if len(grid)-steps<n: raise ValueError("Invalid common lag range.")
            for start in np.arange(0.,t[-1]-win+1e-8,win):
                tt=start+grid
                prog=(np.interp(start+win/2,t,O)-O[0])/span if span>0 else np.nan
                for fname,q in flows.items():
                    qg=np.interp(tt,t,q); x0=qg[:n]
                    if np.std(x0)<1e-8: continue
                    for gas,g in gases.items():
                        y0=np.interp(tt,t,g)
                        if np.std(y0)<1e-8: continue
                        for pw in WHITEN:
                            x,y,phi=prepare(x0,y0,pw)
                            if len(x)<5 or np.std(x)<1e-9 or np.std(y)<1e-9: continue
                            for metric in METRICS:
                                for mode in PEAKS:
                                    lag,r,ar,bound,no_pos=peak(x,y,metric,mode)
                                    rows.append({
                                      "scenario":scen,"heat":heat,"gas":gas,"window_s":win,
                                      "start_s":float(start),"flow_method":fname,
                                      "flow_mean_m3_s":float(np.mean(x0)),"oxygen_progress":prog,
                                      "metric":metric,"prewhiten":pw,"ar1_phi":phi,
                                      "peak_mode":mode,"lag_s":lag,"selected_signed_r":r,
                                      "peak_abs_r":ar,"boundary_peak":bound,
                                      "no_positive_peak":no_pos,"search_max_s":MAX_LAG,
                                    })
    return pd.DataFrame(rows)

def boot_ci(x,seed,n=5000):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    if len(x)<3: return np.nan,np.nan
    rng=np.random.default_rng(seed)
    d=rng.choice(x,size=(n,len(x)),replace=True)
    m=np.median(d,axis=1)
    return float(np.quantile(m,.025)),float(np.quantile(m,.975))

def summarize(w):
    keys=["scenario","gas","window_s","flow_method","metric","prewhiten","peak_mode"]
    ss=[]; cs=[]; hs=[]
    for key,g in w.groupby(keys,dropna=False):
        scen,gas,win,method,metric,pw,mode=key
        v=g[np.isfinite(g.lag_s)]
        ss.append({
          "scenario":scen,"gas":gas,"window_s":win,"flow_method":method,
          "metric":metric,"prewhiten":pw,"peak_mode":mode,
          "n_windows":len(g),"n_valid":len(v),"n_heats":g.heat.nunique(),
          "median_lag_s":v.lag_s.median() if len(v) else np.nan,
          "q25_lag_s":v.lag_s.quantile(.25) if len(v) else np.nan,
          "q75_lag_s":v.lag_s.quantile(.75) if len(v) else np.nan,
          "median_signed_r":v.selected_signed_r.median() if len(v) else np.nan,
          "median_abs_peak_r":v.peak_abs_r.median() if len(v) else np.nan,
          "boundary_rate":v.boundary_peak.mean() if len(v) else np.nan,
          "no_positive_rate":g.no_positive_peak.mean(),
        })
        ds=[]
        for heat,h in v.groupby("heat"):
            if len(h)<4: continue
            q1,q2=h.flow_mean_m3_s.quantile([1/3,2/3]).to_numpy()
            lo=h[h.flow_mean_m3_s<=q1]; hi=h[h.flow_mean_m3_s>q2]
            if len(lo)==0 or len(hi)==0: continue
            d=float(hi.lag_s.median()-lo.lag_s.median()); ds.append(d)
            hs.append({"scenario":scen,"gas":gas,"window_s":win,"flow_method":method,
              "metric":metric,"prewhiten":pw,"peak_mode":mode,"heat":heat,
              "n_low":len(lo),"n_high":len(hi),"low_lag_s":lo.lag_s.median(),
              "high_lag_s":hi.lag_s.median(),"high_minus_low_lag_s":d})
        seed=20261128+int(win)+(0 if scen=="H0" else 100)+(0 if gas=="CO2" else 10)+(0 if method=="analytic" else 1)+(0 if metric=="pearson" else 2)+(0 if pw else 4)+(0 if mode=="absolute" else 8)
        low,high=boot_ci(ds,seed)
        cs.append({"scenario":scen,"gas":gas,"window_s":win,"flow_method":method,
          "metric":metric,"prewhiten":pw,"peak_mode":mode,"n_heats":len(ds),
          "median_high_minus_low_lag_s":float(np.median(ds)) if ds else np.nan,
          "mean_high_minus_low_lag_s":float(np.mean(ds)) if ds else np.nan,
          "n_heats_high_flow_shorter":int(np.sum(np.asarray(ds)<0)),
          "cluster_bootstrap_95_low":low,"cluster_bootstrap_95_high":high})
    return pd.DataFrame(ss),pd.DataFrame(cs),pd.DataFrame(hs)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,default=Path("results/q2/final_freeze_2026-11-28/final_curves.csv"))
    ap.add_argument("--output",type=Path,default=Path("results/q2/lag_round2_2026-11-28"))
    a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    raw=a.input.read_bytes(); d=pd.read_csv(a.input)
    req=["scenario","heat","time_s","oxygen_cumulative_raw_m3","oxygen_cumulative_fit_m3",
         "oxygen_flow_analytic_m3_s","CO2_flow_k2_m3_s","CO_flow_k2_m3_s"]
    if len(d)!=27178 or d[req].isna().any().any(): raise ValueError("Expected 27,178 complete frozen Q2.1 rows.")
    if set(d.scenario)!=set(SCENARIOS) or set(d.heat)!=set(HEATS): raise ValueError("Unexpected heat/scenario set.")
    w=build(d); s,c,h=summarize(w)
    w.to_csv(a.output/"window_level.csv.gz",index=False,compression="gzip")
    s.to_csv(a.output/"config_summary.csv",index=False,float_format="%.8g")
    c.to_csv(a.output/"flow_contrast_summary.csv",index=False,float_format="%.8g")
    h.to_csv(a.output/"heat_flow_contrasts.csv",index=False,float_format="%.8g")
    manifest={"input_sha256":hashlib.sha256(raw).hexdigest(),"rows":len(d),"heats":list(HEATS),
      "scenarios":list(SCENARIOS),"windows_s":list(WINDOWS),"common_lag_search_s":[0,MAX_LAG],
      "grid_s":GRID,"input_segment":"first approximately two-thirds of each window",
      "response_segment":"full fixed window","flow_methods":list(flow_methods(d[(d.scenario=="H0")&(d.heat==HEATS[0])].sort_values("time_s"))),
      "metrics":list(METRICS),"prewhiten":list(WHITEN),"peak_modes":list(PEAKS),
      "independent_unit":"heat","bootstrap_unit":"heat","bootstrap_replicates":5000,
      "primary":{"scenario":"H0","window_s":90,"flow_method":"analytic","metric":"pearson","prewhiten":True,"peak_mode":"absolute"},
      "raw_source_modified":False,"prediction_models":0}
    (a.output/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    primary=s[(s.scenario=="H0")&(s.flow_method=="analytic")&(s.metric=="pearson")&(s.prewhiten)&(s.peak_mode=="absolute")].sort_values(["gas","window_s"])
    contrast=c[(c.scenario=="H0")&(c.flow_method=="analytic")&(c.metric=="pearson")&(c.prewhiten)&(c.peak_mode=="absolute")].sort_values(["gas","window_s"])
    lines=["# Q2.3 第二轮实验：统一搜索范围与时滞口径敏感性","",
      "范围：仅氧枪供氧变化与 CO/CO2 烟气响应有效时滞；不涉及第三问。","",
      f"输入 SHA-256：{manifest['input_sha256']}","",
      "## 设计","",
      "- 90/120/180 s 窗口统一搜索 0–30 s，2 s 网格；氧枪输入取窗口前约 2/3，烟气响应取整窗。",
      "- 比较 Pearson/Spearman、去趋势/AR(1) 输入预白化、正相关峰/绝对相关峰。",
      "- 比较 Q2.1 解析导数与累计耗氧原始/拟合曲线的 2/30/60 s 中心差分。",
      "- 高低流量比较在炉内按流量三分位数分组，Bootstrap 以炉次为单位。","",
      "## 主口径结果（H0、解析导数、Pearson、预白化、绝对峰）","",
      "| 气体 | 窗口 | 有效窗口 | 中位时滞 | IQR | 中位绝对相关峰 | 边界峰比例 |",
      "|---|---:|---:|---:|---:|---:|---:|"]
    for r in primary.itertuples():
        lines.append(f"| {r.gas} | {r.window_s}s | {r.n_valid}/{r.n_windows} | {r.median_lag_s:.1f}s | {r.q25_lag_s:.1f}–{r.q75_lag_s:.1f}s | {r.median_abs_peak_r:.3f} | {r.boundary_rate:.1%} |")
    lines += ["","## 主口径炉内高低流量对照","","| 气体 | 窗口 | 炉数 | 高流量减低流量 | 炉次 Bootstrap 95%区间 | 高流量更短炉数 |",
      "|---|---:|---:|---:|---:|---:|"]
    for r in contrast.itertuples():
        lines.append(f"| {r.gas} | {r.window_s}s | {r.n_heats} | {r.median_high_minus_low_lag_s:+.1f}s | [{r.cluster_bootstrap_95_low:+.1f},{r.cluster_bootstrap_95_high:+.1f}]s | {r.n_heats_high_flow_shorter}/{r.n_heats} |")
    lines += ["","## 解释边界","","- 统一搜索上限排除了原 90/180 s 最大滞后范围不同造成的直接混杂；窗口长度本身仍会改变输入段长度和过程阶段覆盖。",
      "- 正相关峰与绝对峰不同，意味着单一绝对峰可能混合不同方向的响应；必须同时报告峰的符号。",
      "- AR(1) 预白化按输入流量估计滤波系数，并将同一滤波器施加于烟气序列；是交叉相关诊断，不是物理传递函数。",
      "- 有限差分只作流量估计敏感性对照，Q2.1 解析导数仍为主定义。",
      "- 十炉是独立单位；配置扫描是敏感性分析，不按逐窗口数伪造样本量。",
      "- 本轮不拟合新的复杂时滞函数，不预设高流量必须对应更短时滞。","",
      "文件：config_summary.csv、flow_contrast_summary.csv、heat_flow_contrasts.csv、window_level.csv.gz、manifest.json。",""]
    (a.output/"REPORT.md").write_text("\n".join(lines),encoding="utf-8")
    print("Q2.3 round-2 complete; windows:",len(w),"configs:",len(s))
    print(primary.to_string(index=False))
    print(contrast.to_string(index=False))

if __name__=="__main__":
    main()
