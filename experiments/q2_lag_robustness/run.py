#!/usr/bin/env python3
"""Q2.2 effective-lag robustness checks; descriptive only, no prediction."""
from __future__ import annotations
import argparse, hashlib, json
from math import comb
from pathlib import Path
import numpy as np
import pandas as pd

HEATS=["6B03205","6B03206","6B03270","6B03271","6C02057","6C02058","6C02059","6C02060","6C02107","6C02108"]
SCENARIOS=["H0","H1"]
WINDOWS=[90,180]

def cdiff(t, y, half):
    a=np.maximum(t-half,t[0]); b=np.minimum(t+half,t[-1])
    return (np.interp(b,t,y)-np.interp(a,t,y))/(b-a)

def methods(f):
    t=f.time_s.to_numpy(float)
    return {
      "analytic":f.oxygen_flow_analytic_m3_s.to_numpy(float),
      "raw_cum_30s":cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),15.),
      "raw_cum_60s":cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),30.),
      "fit_cum_30s":cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),15.),
      "fit_cum_60s":cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),30.),
    }

def detrend(x):
    x=np.asarray(x,float); t=np.arange(len(x),dtype=float)
    X=np.column_stack([np.ones(len(x)),t])
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

def whiten(x,y):
    x=detrend(x); y=detrend(y)
    den=x[:-1]@x[:-1]
    phi=np.clip((x[:-1]@x[1:])/den if den>0 else 0.,-.98,.98)
    return x[1;]-phi*x[:-1],y[1:]-phi*y[:-1]

def corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x-x.mean(); y=y-y.mean(); den=np.sqrt((x@x)*(y@y))
    return float(x@y/den) if den>0 else np.nan

def bestlag(x,y,maxlag):
    a=[abs(corr(x,y[k:k+len(x)])) for k in range(maxlag+1)]
    k=int(np.nanargmax(a))
    return 2*k,float(a[k])

def make_windows(f,q,y,scen,heat,method,signal,win,offset,kind,true_lag=np.nan,noise=np.nan):
    t=f.time_s.to_numpy(float); O=f.oxygen_cumulative_fit_m3.to_numpy(float)
    rows=[]
    for start in np.arange(offset,t[-1]-win+1e-8,win):
        grid=start+2*np.arange(win//2); qr=np.interp(grid,t,q); yr=np.interp(grid,t,y)
        n=win//3
        if np.std(qr[:n])<1e-7 or np.std(yr)<1e-7: continue
        xw,yw=whiten(qr[:n],yr)
        if np.std(xw)<1e-9 or np.std(yw)<1e-9: continue
        lag,score=bestlag(xw,yw,win//6)
        progress=(np.interp(start+win/2,t,O)-O[0])/(O[-1]-O[0])
        rows.append(dict(kind=kind,scenario=scen,heat=heat,method=method,signal=signal,window_s=win,
            offset_s=offset,start_s=start,flow_m3_s=float(qr[:n].mean()),oxygen_progress=progress,
            estimated_lag_s=lag,peak_abs_correlation=score,true_lag_s=true_lag,noise_fraction=noise))
    return rows

def upper_sign_p(k,n):
    return sum(comb(n,j) for j in range(k,n+1))/(2**n) if n else np.nan

def summarize(w):
    keys=["kind","scenario","method","signal","window_s","offset_s"]
    out=[]
    for key,g in w.groupby(keys,dropna=False):
        kind,scen,method,signal,win,offset=key
        ds=[]
        for heat,h in g.groupby("heat"):
            if len(h)<4: continue
            q1,q2=np.quantile(h.flow_m3_s,[1/3,2/3])
            lo=h[h.flow_m3_s<=q1]; hi=h[h.flow_m3_s>q2]
            if len(lo) and len(hi): ds.append((heat,float(hi.estimated_lag_s.median()-lo.estimated_lag_s.median())))
        d=np.array([x[1] for x in ds],float); nz=d[d!=0]; k=int((nz<0).sum())
        parts=[]
        for heat,h in g.groupby("heat"):
            x=h.flow_m3_s.to_numpy(float); y=h.estimated_lag_s.to_numpy(float); s=h.oxygen_progress.to_numpy(float)
            if len(x)<3: continue
            X=np.column_stack([np.ones(len(s)),s])
            xr=x-X@np.linalg.lstsq(X,x,rcond=None)[0]
            yr=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
            parts.append((heat,float(xr@xr/len(xr)),float(xr@yr/len(yr))))
        den=sum(p[1] for p in parts); num=sum(p[2] for p in parts)
        slope=num/den if den>0 else np.nan
        loo=[]
        for omit,_,_ in parts:
            rest=[p for p in parts if p[0]!=omit]; dd=sum(p[1] for p in rest)
            if dd>0: loo.append(sum(p[2] for p in rest)/dd)
        out.append(dict(kind=kind,scenario=scen,method=method,signal=signal,window_s=win,offset_s=offset,
          n_windows=len(g),n_heats=g.heat.nunique(),median_lag_s=g.estimated_lag_s.median(),
          lag_q25_s=g.estimated_lag_s.quantile(.25),lag_q75_s=g.estimated_lag_s.quantile(.75),
          median_high_minus_low_lag_s=np.median(d) if len(d) else np.nan,n_heats_contrast=len(d),
          n_heats_high_flow_shorter=int((d<0).sum()) if len(d) else np.nan,
          n_nonzero_heat_contrasts=len(nz),
          sign_test_one_sided_p=upper_sign_p(k,len(nz)) if kind=="actual" else np.nan,
          stage_linear_adjusted_slope=slope,slope_loo_min=min(loo) if loo else np.nan,
          slope_loo_max=max(loo) if loo else np.nan,n_loo_negative=sum(x<0 for x in loo),n_loo=len(loo)))
    return pd.DataFrame(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,default=Path("results/q2/final_freeze_2026-11-28/final_curves.csv"))
    ap.add_argument("--output",type=Path,default=Path("results/q2/lag_hypothesis_robustness"))
    a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    raw=a.input.read_bytes(); d=pd.read_csv(a.input)
    required=["scenario","heat","time_s","oxygen_cumulative_raw_m3","oxygen_cumulative_fit_m3",
      "oxygen_flow_analytic_m3_s","CO2_flow_k2_m3_s","CO_flow_k2_m3_s"]
    if len(d)!=27178 or d[required].isna().any().any(): raise ValueError("Expected 27,178 complete frozen Q2.1 rows.")
    if set(d.scenario)!=set(SCENARIOS) or set(d.heat)!=set(HEATS): raise ValueError("Unexpected heat/scenario set.")
    actual=[]; synthetic=[]
    for (scen,heat),f in d.groupby(["scenario","heat"]):
        f=f.sort_values("time_s"); t=f.time_s.to_numpy(float); m=methods(f); qtrue=m["analytic"]
        for method,q in m.items():
            for win in WINDOWS:
                for off in [0,win//2]:
                    for gas in ["CO","CO2"]:
                        actual += make_windows(f,q,f[f"{gas}_flow_k2_m3_s"].to_numpy(float),scen,heat,method,gas,win,off,"actual")
        # Deterministic diagnostic: fixed 20-s delay from the analytic flow signal; no added noise.
        tau=20.
        y=np.interp(t-tau,t,qtrue,left=qtrue[0],right=qtrue[-1])
        noise_fraction=0.0  # deterministic injection: isolate estimator/window effects
        for method,q in m.items():
            for win in WINDOWS:
                for off in [0,win//2]:
                    synthetic += make_windows(f,q,y,scen,heat,method,"synthetic_fixed_20s",win,off,"synthetic",tau,noise_fraction)
    aw=pd.DataFrame(actual); sw=pd.DataFrame(synthetic)
    ass=summarize(aw); sss=summarize(sw)
    # Family-wise correction for the 40 primary tests: 5 flow estimators x
    # 2 scenarios x 2 gases x 2 window lengths, using unshifted windows.
    def holm_adjust(values):
        values=np.asarray(values,float); order=np.argsort(values); adjusted=np.empty(len(values)); running=0.0
        for rank,index in enumerate(order):
            running=max(running,(len(values)-rank)*values[index])
            adjusted[index]=min(1.0,running)
        return adjusted
    primary=ass.offset_s==0
    ass["sign_test_holm_p_primary"]=np.nan
    ass.loc[primary,"sign_test_holm_p_primary"]=holm_adjust(
        ass.loc[primary,"sign_test_one_sided_p"].fillna(1.0).to_numpy()
    )
    aw.to_csv(a.output/"actual_lag_windows.csv",index=False,float_format="%.8g")
    ass.to_csv(a.output/"actual_lag_summary.csv",index=False,float_format="%.8g")
    sw.to_csv(a.output/"synthetic_fixed_lag_windows.csv",index=False,float_format="%.8g")
    sss.to_csv(a.output/"synthetic_fixed_lag_summary.csv",index=False,float_format="%.8g")
    manifest={"input_sha256":hashlib.sha256(raw).hexdigest(),"input_rows":len(d),"heats":HEATS,"scenarios":SCENARIOS,
      "flow_methods":list(methods(d[(d.scenario=="H0")&(d.heat==HEATS[0])].sort_values("time_s"))),
      "windows_s":WINDOWS,"offsets_s":{"90":[0,45],"180":[0,90]},"grid_s":2,"prediction_models":0,
      "synthetic_fixed_lag_test":{"response":"Q2.1 analytic oxygen-flow signal delayed by exactly 20 s",
        "candidate_inputs":"five oxygen-flow estimators","added_noise_sd":"0 (deterministic injection)",
        "interpretation":"estimator/window diagnostic only; not a physical gas-process simulation"},
      "raw_source_modified":False}
    (a.output/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("Wrote Q2.2 lag robustness results:",a.output)
    print(ass[(ass.method=="analytic")&(ass.window_s==90)].to_string(index=False))
    print(sss[(sss.window_s==90)&(sss.offset_s==0)][["scenario","method","median_lag_s","stage_linear_adjusted_slope","slope_loo_min","slope_loo_max"]].to_string(index=False))
if __name__=="__main__": main()
