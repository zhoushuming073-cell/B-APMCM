#!/usr/bin/env python3
"""Descriptive relationship and effective-lag study for APMCM Problem B Q2.
No prediction models are trained or evaluated. Input is the frozen Q2.1 final_curves.csv.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

HEATS=["6B03205","6B03206","6B03270","6B03271","6C02057","6C02058","6C02059","6C02060","6C02107","6C02108"]
GASES=["CO2","CO"]
SCENARIOS=["H0","H1"]
WINDOWS=[90,180]
METRICS=["pearson","spearman"]

def corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x-x.mean(); y=y-y.mean(); den=np.sqrt(x@x*(y@y))
    return float(x@y/den) if den>0 else np.nan

def assoc(x,y,metric):
    if metric=="pearson": return corr(x,y)
    return corr(pd.Series(x).rank(method="average"),pd.Series(y).rank(method="average"))

def detrend(x):
    x=np.asarray(x,float); t=np.arange(len(x),dtype=float)
    X=np.column_stack([np.ones(len(x)),t])
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

def whiten(x,y):
    x=detrend(x); y=detrend(y)
    den=x[:-1]@x[:-1]
    phi=np.clip((x[:-1]@x[1:])/den if den>0 else 0.,-.98,.98)
    return x[1:]-phi*x[:-1],y[1:]-phi*y[:-1],float(phi)

def best_lag(x,y,maxlag,metric):
    vals=[abs(assoc(x,y[k:k+len(x)],metric)) for k in range(maxlag+1)]
    k=int(np.nanargmax(vals))
    return k,float(vals[k]),float(assoc(x,y[k:k+len(x)],metric))

def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); z=np.empty(n); run=1.
    for j in range(n-1,-1,-1):
        run=min(run,p[order[j]]*n/(j+1)); z[j]=min(1.,run)
    q=np.empty(n); q[order]=z
    return q

def lag_tables(df):
    rows=[]
    for metric in METRICS:
      for scen in SCENARIOS:
       for win in WINDOWS:
        nin=win//3; maxlag=win//6
        for heat in HEATS:
         f=df[(df.scenario==scen)&(df.heat==heat)].sort_values("time_s")
         t=f.time_s.to_numpy(float); q=f.oxygen_flow_analytic_m3_s.to_numpy(float)
         O=f.oxygen_cumulative_fit_m3.to_numpy(float)
         for start in np.arange(0,t[-1]-win+1e-8,win):
          grid=start+2*np.arange(win//2); xr=np.interp(grid[:nin],t,q)
          progress=(np.interp(start+win/2,t,O)-O[0])/(O[-1]-O[0])
          for gas in GASES:
           yr=np.interp(grid,t,f[f"{gas}_flow_k2_m3_s"].to_numpy(float))
           if np.std(xr)<1e-7 or np.std(yr)<1e-7: continue
           x,y,phi=whiten(xr,yr)
           if np.std(x)<1e-9 or np.std(y)<1e-9: continue
           k,peak,r=best_lag(x,y,maxlag,metric); null=[]
           for shift in range(len(y)):
            if shift<=maxlag or shift>=len(y)-maxlag: continue
            null.append(best_lag(x,np.roll(y,-shift),maxlag,metric)[1])
           p=(1+sum(v>=peak for v in null))/(1+len(null))
           rows.append(dict(metric=metric,scenario=scen,window_s=win,heat=heat,gas=gas,start_s=start,
             flow_m3_s=xr.mean(),oxygen_progress=progress,ar1_phi=phi,lag_s=2*k,peak_score=peak,
             signed_association=r,p_value=p,n_null=len(null),boundary_peak=k==maxlag))
    w=pd.DataFrame(rows); w["bh_q"]=np.nan
    for _,ix in w.groupby(["metric","scenario","window_s","gas"]).groups.items():
        w.loc[ix,"bh_q"]=bh(w.loc[ix,"p_value"].to_numpy())
    sums=[]; bins=[]; slopes=[]; contrasts=[]; contrast_sums=[]
    for (metric,scen,win,gas),g in w.groupby(["metric","scenario","window_s","gas"]):
      sums.append(dict(metric=metric,scenario=scen,window_s=win,gas=gas,n_windows=len(g),n_heats=g.heat.nunique(),
        median_lag_s=g.lag_s.median(),lag_q25_s=g.lag_s.quantile(.25),lag_q75_s=g.lag_s.quantile(.75),
        median_peak_score=g.peak_score.median(),supported_windows_q_le_0_1=int((g.bh_q<=.1).sum()),
        boundary_rate=g.boundary_peak.mean()))
      v=g.flow_m3_s.to_numpy(float); c=np.quantile(v,[1/3,2/3])
      labels=np.where(v<=c[0],"low",np.where(v<=c[1],"middle","high"))
      z=g.copy(); z["flow_bin"]=labels
      for label,b in z.groupby("flow_bin"):
        bins.append(dict(metric=metric,scenario=scen,window_s=win,gas=gas,flow_bin=label,n_windows=len(b),n_heats=b.heat.nunique(),
          flow_min=b.flow_m3_s.min(),flow_median=b.flow_m3_s.median(),flow_max=b.flow_m3_s.max(),
          lag_median_s=b.lag_s.median(),lag_q25_s=b.lag_s.quantile(.25),lag_q75_s=b.lag_s.quantile(.75),
          boundary_rate=b.boundary_peak.mean(),supported_windows=int((b.bh_q<=.1).sum())))
      for degree,name in [(0,"within_heat"),(1,"stage_linear_adjusted"),(2,"stage_quadratic_adjusted")]:
        parts=[]
        for heat,h in g.groupby("heat"):
          x=h.flow_m3_s.to_numpy(float); y=h.lag_s.to_numpy(float); s=h.oxygen_progress.to_numpy(float)
          if len(x)<degree+2: continue
          if degree:
            X=np.column_stack([s**k for k in range(degree+1)])
            x=x-X@np.linalg.lstsq(X,x,rcond=None)[0]
            y=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
          parts.append((x@x/len(x),x@y/len(x)))
        xx=sum(a for a,b in parts); xy=sum(b for a,b in parts)
        slope=xy/xx if xx>0 else np.nan
        rng=np.random.default_rng(20261128+degree+win+(0 if scen=="H0" else 100)+(0 if gas=="CO2" else 10)+(0 if metric=="pearson" else 1))
        boot=[]; a=np.asarray(parts)
        if len(a)>=3:
          for _ in range(3000):
            zz=a[rng.integers(0,len(a),len(a))].sum(axis=0)
            if zz[0]>0: boot.append(zz[1]/zz[0])
        ci=np.quantile(boot,[.025,.975]) if boot else [np.nan,np.nan]
        slopes.append(dict(metric=metric,scenario=scen,window_s=win,gas=gas,adjustment=name,
          slope_s_per_m3s=slope,cluster_bootstrap_95_low=ci[0],cluster_bootstrap_95_high=ci[1],
          n_windows=len(g),n_heats=g.heat.nunique(),n_supported=int((g.bh_q<=.1).sum())))
      ds=[]
      for heat,h in g.groupby("heat"):
        if len(h)<4: continue
        q1,q2=np.quantile(h.flow_m3_s,[1/3,2/3]); lo=h[h.flow_m3_s<=q1]; hi=h[h.flow_m3_s>q2]
        if lo.empty or hi.empty: continue
        d=float(hi.lag_s.median()-lo.lag_s.median()); ds.append(d)
        contrasts.append(dict(metric=metric,scenario=scen,window_s=win,gas=gas,heat=heat,
          q_low_median=lo.flow_m3_s.median(),lag_low_median=lo.lag_s.median(),
          q_high_median=hi.flow_m3_s.median(),lag_high_median=hi.lag_s.median(),
          high_minus_low_lag_s=d,n_low=len(lo),n_high=len(hi)))
      if ds:
        a=np.asarray(ds); rng=np.random.default_rng(20261128+win+(0 if scen=="H0" else 100)+(0 if gas=="CO2" else 10)+(0 if metric=="pearson" else 1))
        boot=np.asarray([np.median(a[rng.integers(0,len(a),len(a))]) for _ in range(5000)])
        contrast_sums.append(dict(metric=metric,scenario=scen,window_s=win,gas=gas,n_heats=len(a),
          median_high_minus_low_lag_s=np.median(a),mean_high_minus_low_lag_s=np.mean(a),
          heats_with_shorter_lag_at_high_flow=int((a<0).sum()),
          cluster_bootstrap_95_low=np.quantile(boot,.025),cluster_bootstrap_95_high=np.quantile(boot,.975)))
    return w,pd.DataFrame(sums),pd.DataFrame(bins),pd.DataFrame(slopes),pd.DataFrame(contrasts),pd.DataFrame(contrast_sums)

def relationship_tables(df):
    dose_rows=[]; stage=[]; yields=[]; blocks=[]; totals=[]
    for scen in SCENARIOS:
      frames={h:df[(df.scenario==scen)&(df.heat==h)].sort_values("time_s") for h in HEATS}
      xmax=min(f.oxygen_cumulative_fit_m3.iloc[-1]-f.oxygen_cumulative_fit_m3.iloc[0] for f in frames.values())
      grid=np.linspace(0,xmax,101); edges=np.linspace(0,xmax,5)
      for heat,f in frames.items():
        t=f.time_s.to_numpy(float); O=f.oxygen_cumulative_fit_m3.to_numpy(float); x=O-O[0]
        keep=np.r_[True,np.diff(x)>1e-10]; xu=x[keep]; tu=t[keep]; tg=np.interp(grid,xu,tu)
        te=np.interp(edges,xu,tu)
        for gas in GASES:
          flow=f[f"{gas}_flow_k2_m3_s"].to_numpy(float); E=f[f"{gas}_cumulative_k2_m3"].to_numpy(float)
          fg=np.interp(tg,t,flow); eg=np.interp(tg,t,E)
          for dose_val,tm,g,e in zip(grid,tg,fg,eg):
            dose_rows.append(dict(scenario=scen,heat=heat,gas=gas,oxygen_dose_m3=dose_val,
              oxygen_progress=dose_val/xmax,time_s=tm,gas_flow_m3_s=g,cumulative_emission_m3=e))
          ee=np.interp(te,t,E)
          for j in range(4):
            m=(t>=te[j])&(t<=te[j+1]); avg=float(flow[m].mean()) if m.any() else np.nan
            de=float(ee[j+1]-ee[j]); do=float(edges[j+1]-edges[j])
            stage.append(dict(scenario=scen,heat=heat,gas=gas,stage=j+1,mean_gas_flow_m3_s=avg))
            yields.append(dict(scenario=scen,heat=heat,gas=gas,stage=j+1,delta_oxygen_m3=do,
              delta_emission_m3=de,marginal_emission_per_oxygen=de/do))
          totals.append(dict(scenario=scen,heat=heat,gas=gas,net_oxygen_m3=x[-1],total_emission_m3=E[-1],
            emission_per_oxygen=E[-1]/x[-1],mean_gas_flow_m3_s=flow.mean()))
          do=[]; gf=[]; de=[]
          for a in np.arange(0,t[-1]-30+1e-8,30.):
            b=a+30; oab=np.interp([a,b],t,O); eab=np.interp([a,b],t,E); m=(t>=a)&(t<=b)
            if m.sum()>1:
              do.append(oab[1]-oab[0]); gf.append(flow[m].mean()); de.append(eab[1]-eab[0])
          blocks.append(dict(scenario=scen,heat=heat,gas=gas,n_blocks=len(do),
            spearman_dO_vs_mean_gas_flow=assoc(do,gf,"spearman"),spearman_dO_vs_dE=assoc(do,de,"spearman")))
    dose=pd.DataFrame(dose_rows); ds=[]
    for (scen,gas),g in dose.groupby(["scenario","gas"]):
      for d,b in g.groupby("oxygen_dose_m3"):
        e=b.cumulative_emission_m3.to_numpy(); f=b.gas_flow_m3_s.to_numpy()
        ds.append(dict(scenario=scen,gas=gas,oxygen_dose_m3=d,emission_median_m3=np.median(e),
          emission_q25_m3=np.quantile(e,.25),emission_q75_m3=np.quantile(e,.75),
          gas_flow_median_m3_s=np.median(f),gas_flow_q25_m3_s=np.quantile(f,.25),
          gas_flow_q75_m3_s=np.quantile(f,.75),n_heats=len(b)))
    return dose,pd.DataFrame(ds),pd.DataFrame(stage),pd.DataFrame(yields),pd.DataFrame(blocks),pd.DataFrame(totals)

def synthetic_check():
    x=np.sin(np.linspace(0,12*np.pi,180))+.3*np.sin(np.linspace(0,31*np.pi,180)); out=[]
    for delay in [0,6,12,20,30]:
      k=delay//2; y=np.r_[np.repeat(x[0],k),x[:-k]] if k else x.copy()
      a,b,_=whiten(x[:90],y[:120]); lag,score,_=best_lag(a,b,15,"pearson")
      out.append(dict(injected_delay_s=delay,recovered_delay_s=2*lag,peak_abs_r=score,recovered_exactly=2*lag==delay))
    return pd.DataFrame(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=Path('results/q2/final_freeze_2026-11-28/final_curves.csv'))
    ap.add_argument('--output',type=Path,default=Path('results/q2/relationship_study_final'))
    a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=True); d=pd.read_csv(a.input)
    required=['scenario','heat','time_s','oxygen_cumulative_fit_m3','oxygen_flow_analytic_m3_s',
      'CO2_flow_k2_m3_s','CO_flow_k2_m3_s','CO2_cumulative_k2_m3','CO_cumulative_k2_m3']
    if len(d)!=27178 or d[required].isna().any().any() or set(d.scenario)!={'H0','H1'} or set(d.heat)!=set(HEATS):
      raise ValueError('Expected 27,178 complete rows from frozen Q2.1 curves (10 heats x H0/H1).')
    for (s,h),g in d.groupby(['scenario','heat']):
      if not np.all(np.diff(g.sort_values('time_s').time_s.to_numpy())>0): raise ValueError(f'Invalid time axis: {s}/{h}')
    w,ls,fb,sl,wh,whs=lag_tables(d); dc,ds,st,yi,bl,tot=relationship_tables(d); sc=synthetic_check()
    out={'lag_windows.csv':w,'lag_summary.csv':ls,'lag_flow_bins.csv':fb,'lag_flow_slopes.csv':sl,
      'within_heat_flow_contrasts.csv':wh,'within_heat_flow_contrast_summary.csv':whs,'dose_response_curves.csv':dc,
      'dose_response_summary.csv':ds,'gas_flow_by_oxygen_stage.csv':st,'marginal_emission_yield_by_stage.csv':yi,
      'block_relationships.csv':bl,'heat_emission_totals.csv':tot,'synthetic_lag_recovery.csv':sc}
    for name,frame in out.items(): frame.to_csv(a.output/name,index=False,float_format='%.10g')
    manifest={'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'rows':len(d),'heats':10,
      'scenarios':SCENARIOS,'oxygen_signal':'Q2.1 nonnegative integral B-spline analytic derivative',
      'gas_signal':'Q2.1 k2/PCHIP flow and exact cumulative integral','prediction_models':0,
      'primary_window_s':90,'sensitivity_window_s':180,'lag_metrics':METRICS,
      'stage_adjustment':['none','linear','quadratic'],'raw_source_modified':False}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Q2 final descriptive relationship study complete; prediction models: 0')
    print(ls[ls.metric=='pearson'].to_string(index=False))
    print(whs[(whs.metric=='pearson')&(whs.window_s==90)].to_string(index=False))

if __name__=='__main__': main()
