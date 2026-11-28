from __future__ import annotations
import argparse, hashlib, json, math, platform, sys, time
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linear_sum_assignment
from scipy.stats import t as student_t, binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HEATS=['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108']
GASES=['CO','CO2']; SCENARIOS=['H0','H1']; WINDOWS=[90,180]
SEED=20261128; N_SIM=100; DT_H0=.5

CONFIG={
 'main_scenario':'H0','main_window_s':90,'sensitivity_window_s':180,'window_step_s':'equal to window length (non-overlapping)',
 'oxygen_input_segment':'first one-third of window duration (60 s for 90 s; 120 s for 180 s), sampled every 2 s',
 'gas_observation_segment':'full window', 'sampling_grid_s':2,'lag_grid_s':2,
 'max_lag_s':{'90s':30,'180s':60},'lag_sign':'positive lag means oxygen input precedes gas observation',
 'detrending':'OLS linear trend removed separately from oxygen and gas vectors',
 'prewhitening':'AR(1) phi estimated from detrended oxygen input and applied to both oxygen and gas; phi clipped to [-0.98,0.98]',
 'peak_rule':'maximum absolute Pearson correlation over nonnegative lags; first/smallest lag wins exact ties',
 'null_test':'circularly shifted gas whitened vector; allowed shifts exclude <= maxlag and >= len(y)-maxlag; p=(1+exceedances)/(1+n_null); BH within metric/scenario/window/gas',
 'flow':'Q2.1 nonnegative integral cubic B-spline analytic derivative; mean over same first-2/3 input interval',
 'gas':'Q2.1 k2 local 2-s window means at 2-s knots, PCHIP reconstruction; raw units m3/s per official problem statement',
 'progress_old':'cumulative oxygen at full-window midpoint, normalized by heat total',
 'progress_new':'mean normalized cumulative oxygen over exactly the first-2/3 oxygen-input grid points',
 'simulation':'raw gas observation layer; y_syn(t)=T_gas(t)+beta*sign*Q_fluct(t-tau)+epsilon_block(t), beta=sign*SNR*sd(epsilon)/sd(Q_fluct over valid t>=tau); then same k2/PCHIP and lag estimator',
 'trend_smoother':'centered 60-s moving average (121 rows at H0 0.5-s sampling); same smoother applied to oxygen flow and gas signal',
 'noise_blocks':'other-heats same-gas detrended raw residuals; moving blocks; block length = median first 1/e ACF crossing across heats, rounded up; each donor block centered and final sequence rescaled to target-heat residual SD',
 'snr_levels':[1.0,0.5],'injected_lags_s':[0,10,20,30],'directions':['positive','negative'],
 'replicates_per_random_condition':100,'random_seed':SEED,
 'zero_response_control':'stage trend + empirical residual noise, response coefficient exactly zero; 100 replicates; BH false-positive fraction reported',
 'upper_bound_sensitivity':'same injected delays rerun with 180-s window and 60-s max lag; both window length and max search bound change',
}

def sha256(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); x=x-x.mean(); y=y-y.mean(); den=np.sqrt((x@x)*(y@y))
    return float(x@y/den) if den>0 else np.nan

def assoc(x,y,metric='pearson'):
    if metric=='pearson': return corr(x,y)
    return corr(pd.Series(x).rank(method='average'),pd.Series(y).rank(method='average'))

def detrend(x):
    x=np.asarray(x,float); z=np.arange(len(x),dtype=float); X=np.column_stack([np.ones(len(x)),z])
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

def whiten(x,y):
    x=detrend(x); y=detrend(y); den=x[:-1]@x[:-1]
    phi=np.clip((x[:-1]@x[1:])/den if den>0 else 0.,-.98,.98)
    return x[1:]-phi*x[:-1],y[1:]-phi*y[:-1],float(phi)

def best_lag(x,y,maxlag,metric='pearson'):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if metric=='pearson':
        wins=np.lib.stride_tricks.sliding_window_view(y,len(x))[:maxlag+1]
        xc=x-x.mean();yc=wins-wins.mean(axis=1,keepdims=True)
        den=np.sqrt((xc@xc)*np.sum(yc*yc,axis=1))
        vals=np.divide(yc@xc,den,out=np.full(len(yc),np.nan),where=den>0)
    else:
        vals=np.array([assoc(x,y[k:k+len(x)],metric) for k in range(maxlag+1)],float)
    if not np.isfinite(vals).any(): return None,np.nan,np.nan
    k=int(np.nanargmax(np.abs(vals))); return k,float(abs(vals[k])),float(vals[k])

def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); z=np.empty(n); run=1.
    for j in range(n-1,-1,-1): run=min(run,p[order[j]]*n/(j+1)); z[j]=min(1.,run)
    q=np.empty(n); q[order]=z; return q

def k2_fit(t,y):
    t=np.asarray(t,float); y=np.asarray(y,float)
    knots=np.unique(np.r_[t[0],np.arange(t[0]+2.,t[-1],2.),t[-1]])
    centers=[t[0]]; means=[y[0]]
    for knot in knots[1:-1]:
        take=np.abs(t-knot)<=1.
        if take.sum()>=2: centers.append(float(t[take].mean())); means.append(float(y[take].mean()))
    centers.append(t[-1]); means.append(y[-1])
    return PchipInterpolator(centers,np.maximum(means,0),extrapolate=False)

def get_windows(curves,scenario='H0',window=90,metric='pearson',progress_def='old',compute_p=True,gas_override=None):
    rows=[]; nin=window//3; maxlag=window//6
    for heat in HEATS:
        f=curves[(curves.scenario==scenario)&(curves.heat==heat)].sort_values('time_s')
        if f.empty: raise ValueError(f'missing {scenario}/{heat}')
        t=f.time_s.to_numpy(float); q=f.oxygen_flow_analytic_m3_s.to_numpy(float); O=f.oxygen_cumulative_fit_m3.to_numpy(float)
        total=O[-1]-O[0]
        for start in np.arange(0,t[-1]-window+1e-8,window):
            grid=start+2*np.arange(window//2)
            qx=np.interp(grid[:nin],t,q)
            if progress_def=='old': p=(np.interp(start+window/2,t,O)-O[0])/total
            elif progress_def=='new': p=float(np.mean((np.interp(grid[:nin],t,O)-O[0])/total))
            else: raise ValueError(progress_def)
            for gas in GASES:
                col=f'{gas}_flow_k2_m3_s'
                gy=np.interp(grid,t,(gas_override[(scenario,heat)] if gas_override and (scenario,heat) in gas_override else f[col].to_numpy(float)))
                if np.std(qx)<1e-7 or np.std(gy)<1e-7: continue
                x,y,phi=whiten(qx,gy)
                if np.std(x)<1e-9 or np.std(y)<1e-9: continue
                k,peak,signed=best_lag(x,y,maxlag,metric)
                if k is None: continue
                pval=np.nan; nnull=0
                if compute_p:
                    null=[]
                    for shift in range(len(y)):
                        if shift<=maxlag or shift>=len(y)-maxlag: continue
                        _,sc,_=best_lag(x,np.roll(y,-shift),maxlag,metric)
                        if np.isfinite(sc): null.append(sc)
                    nnull=len(null); pval=(1+sum(v>=peak for v in null))/(1+len(null)) if null else np.nan
                rows.append(dict(metric=metric,scenario=scenario,window_s=window,heat=heat,gas=gas,start_s=float(start),
                    flow_m3_s=float(qx.mean()),oxygen_progress=float(p),lag_s=float(2*k),peak_abs_assoc=peak,
                    signed_assoc_at_peak=signed,ar1_phi=phi,p_value=pval,n_null=nnull,boundary_peak=bool(k==maxlag),
                    zero_lag_peak=bool(k==0),peak_sign='positive' if signed>0 else ('negative' if signed<0 else 'zero')))
    out=pd.DataFrame(rows)
    out['bh_q']=np.nan
    if compute_p and len(out):
        for _,ix in out.groupby(['metric','scenario','window_s','gas']).groups.items(): out.loc[ix,'bh_q']=bh(out.loc[ix,'p_value'].to_numpy(float))
    return out

def flow_bins_table(w):
    rows=[]
    for (metric,scenario,window,gas),g in w.groupby(['metric','scenario','window_s','gas']):
        cut1,cut2=np.quantile(g.flow_m3_s,[1/3,2/3])
        lab=np.where(g.flow_m3_s<=cut1,'low',np.where(g.flow_m3_s<=cut2,'middle','high'))
        z=g.copy();z['flow_bin']=lab
        for name,b in z.groupby('flow_bin'):
            rows.append(dict(metric=metric,scenario=scenario,window_s=window,gas=gas,flow_bin=name,
                cut_low_m3_s=cut1,cut_high_m3_s=cut2,flow_min_m3_s=b.flow_m3_s.min(),flow_max_m3_s=b.flow_m3_s.max(),
                flow_median_m3_s=b.flow_m3_s.median(),n_windows=len(b),n_heats=b.heat.nunique(),
                lag_median_s=b.lag_s.median(),lag_q25_s=b.lag_s.quantile(.25),lag_q75_s=b.lag_s.quantile(.75),
                positive_peak_n=int((b.signed_assoc_at_peak>0).sum()),negative_peak_n=int((b.signed_assoc_at_peak<0).sum()),
                zero_lag_n=int(b.zero_lag_peak.sum()),upper_bound_n=int(b.boundary_peak.sum()),
                bh_q_le_0_1_n=int((b.bh_q<=.1).sum()),bh_denominator=len(b)))
    return pd.DataFrame(rows)

def actual_summary(w):
    rows=[]
    for (metric,scenario,window,gas),g in w.groupby(['metric','scenario','window_s','gas']):
        rows.append(dict(metric=metric,scenario=scenario,window_s=window,gas=gas,n_windows=len(g),n_heats=g.heat.nunique(),
            lag_median_s=g.lag_s.median(),lag_q25_s=g.lag_s.quantile(.25),lag_q75_s=g.lag_s.quantile(.75),
            positive_peak_n=int((g.signed_assoc_at_peak>0).sum()),negative_peak_n=int((g.signed_assoc_at_peak<0).sum()),
            zero_lag_n=int(g.zero_lag_peak.sum()),upper_bound_n=int(g.boundary_peak.sum()),
            bh_q_le_0_1_n=int((g.bh_q<=.1).sum()),bh_denominator=len(g),
            min_empirical_p=float(g.p_value.min()),n_null_min=int(g.n_null.min()),n_null_max=int(g.n_null.max())))
    return pd.DataFrame(rows)

def stage_matrix(p,degree):
    p=np.asarray(p,float);return np.column_stack([p**k for k in range(degree+1)])

def residualize(y,X):
    b=np.linalg.lstsq(X,np.asarray(y,float),rcond=None)[0];return np.asarray(y,float)-X@b,b

def pooled_cluster_reg(g,degree,boundary_mode):
    d=g.copy()
    if boundary_mode=='interior':d=d[~d.boundary_peak].copy()
    d=d.sort_values(['heat','start_s']).reset_index(drop=True);n=len(d); heats=list(d.heat.drop_duplicates())
    if n<20 or len(heats)<3:return None
    stage_cols=[];xres=np.zeros(n)
    for heat in heats:
        ix=np.flatnonzero(d.heat.to_numpy()==heat);h=d.iloc[ix];p=h.oxygen_progress.to_numpy(float);Xp=stage_matrix(p,degree)
        xr,_=residualize(h.flow_m3_s.to_numpy(float),Xp);sd=float(np.std(xr,ddof=1))
        if sd>1e-10:xres[ix]=xr/sd
        for k in range(degree+1):
            col=np.zeros(n);col[ix]=p**k;stage_cols.append(col)
    X=np.column_stack(stage_cols+[xres]);y=d.lag_s.to_numpy(float);xtx=np.linalg.pinv(X.T@X);beta=xtx@X.T@y;res=y-X@beta
    meat=np.zeros((X.shape[1],X.shape[1]))
    groups=d.heat.to_numpy()
    for heat in heats:
        ix=np.flatnonzero(groups==heat);score=X[ix].T@res[ix];meat+=np.outer(score,score)
    G=len(heats);K=X.shape[1];corr=(G/(G-1))*((n-1)/max(n-K,1));cov=corr*xtx@meat@xtx
    se=float(np.sqrt(max(cov[-1,-1],0)));b=float(beta[-1]);crit=float(student_t.ppf(.975,G-1))
    return dict(degree=degree,boundary_mode=boundary_mode,n_windows=n,n_heats=G,beta_s_per_within_heat_sd=b,cluster_se=se,
        ci_low=b-crit*se,ci_high=b+crit*se,df=G-1)

def heat_slopes(g,degree):
    rows=[]
    for heat,h0 in g.groupby('heat'):
        h=h0.sort_values('start_s');p=h.oxygen_progress.to_numpy(float);Xp=stage_matrix(p,degree)
        if len(h)<=degree+2:continue
        x,_=residualize(h.flow_m3_s.to_numpy(float),Xp);sx=float(np.std(x,ddof=1))
        if sx<1e-10:continue
        x=x/sx;y,_=residualize(h.lag_s.to_numpy(float),Xp);slope=float((x@y)/(x@x))
        rows.append(dict(heat=heat,n_windows=len(h),slope_s_per_within_heat_sd=slope))
    return pd.DataFrame(rows)

def match_pairs(g,caliper,progress_col='oxygen_progress'):
    allp=[]
    for heat,h0 in g.groupby('heat',sort=False):
        h=h0.sort_values('start_s').reset_index(drop=True);ql,qh=np.quantile(h.flow_m3_s,[.25,.75])
        lo=h[h.flow_m3_s<ql].sort_values(progress_col).reset_index();hi=h[h.flow_m3_s>qh].sort_values(progress_col).reset_index()
        if lo.empty or hi.empty:continue
        dist=np.abs(hi[progress_col].to_numpy()[:,None]-lo[progress_col].to_numpy()[None,:]);nh,nl=dist.shape
        cost=np.full((nh,nl+nh),caliper+1e-8);cost[:,:nl]=np.where(dist<=caliper,dist,caliper+1.)
        for i in range(nh):cost[i,nl+i]=caliper+1e-8
        rr,cc=linear_sum_assignment(cost)
        for i,j in zip(rr,cc):
            if j>=nl or dist[i,j]>caliper:continue
            H=hi.iloc[i];L=lo.iloc[j]
            allp.append(dict(heat=heat,gas=H.gas,scenario=H.scenario,window_s=H.window_s,metric=H.metric,caliper=caliper,
                high_start_s=H.start_s,low_start_s=L.start_s,high_flow_m3_s=H.flow_m3_s,low_flow_m3_s=L.flow_m3_s,
                flow_diff_m3_s=H.flow_m3_s-L.flow_m3_s,high_progress=H[progress_col],low_progress=L[progress_col],
                progress_gap=abs(H[progress_col]-L[progress_col]),high_lag_s=H.lag_s,low_lag_s=L.lag_s,
                lag_diff_high_minus_low_s=H.lag_s-L.lag_s,high_boundary=H.boundary_peak,low_boundary=L.boundary_peak))
    return pd.DataFrame(allp)

def overlap_diagnostics(g):
    rows=[]
    for heat,h in g.groupby('heat'):
        ql,qh=np.quantile(h.flow_m3_s,[.25,.75]);lo=h[h.flow_m3_s<ql];hi=h[h.flow_m3_s>qh]
        if len(lo) and len(hi):
            pl=(lo.oxygen_progress.min(),lo.oxygen_progress.max());ph=(hi.oxygen_progress.min(),hi.oxygen_progress.max())
            low=max(pl[0],ph[0]);high=min(pl[1],ph[1]);overlap=max(0.,high-low)
            lo_in=int(((lo.oxygen_progress>=low)&(lo.oxygen_progress<=high)).sum()) if overlap>0 else 0
            hi_in=int(((hi.oxygen_progress>=low)&(hi.oxygen_progress<=high)).sum()) if overlap>0 else 0
        else:low=high=overlap=np.nan;lo_in=hi_in=0
        rows.append(dict(heat=heat,n_windows=len(h),q25_flow_m3_s=ql,q75_flow_m3_s=qh,
            low_windows=len(lo),high_windows=len(hi),low_progress_min=pl[0] if len(lo) else np.nan,
            low_progress_max=pl[1] if len(lo) else np.nan,high_progress_min=ph[0] if len(hi) else np.nan,
            high_progress_max=ph[1] if len(hi) else np.nan,overlap_progress_low=low,overlap_progress_high=high,
            overlap_progress_width=overlap,low_windows_inside_overlap=lo_in,high_windows_inside_overlap=hi_in,
            any_progress_overlap=bool(overlap>0) if np.isfinite(overlap) else False))
    return pd.DataFrame(rows)

def rolling_trend(x,window=121):
    return pd.Series(np.asarray(x,float)).rolling(window,center=True,min_periods=1).mean().to_numpy(float)

def acf_efold_seconds(x,dt=.5,max_lag_s=120):
    x=np.asarray(x,float);lim=min(len(x)//2,int(max_lag_s/dt))
    if np.std(x)<1e-12:return np.nan
    for k in range(1,lim):
        a=corr(x[:-k],x[k:])
        if np.isfinite(a) and a<=np.exp(-1):return k*dt
    return lim*dt

def acf1(x):
    x=np.asarray(x,float)
    return corr(x[:-1],x[1:]) if len(x)>2 else np.nan

def moving_block_noise(n,target_heat,residuals,block_len,sigma_target,rng):
    donors=[h for h in residuals if h!=target_heat and len(residuals[h])>block_len+1]
    if not donors: return np.zeros(n)
    out=[]
    while len(out)<n:
        donor=donors[int(rng.integers(0,len(donors)))];r=residuals[donor]
        start=int(rng.integers(0,len(r)-block_len+1));b=r[start:start+block_len].copy()
        out.extend(b.tolist())
    z=np.asarray(out[:n],float);z=z-z.mean();sd=np.std(z)
    return z*(sigma_target/sd) if sd>1e-12 else np.zeros(n)

def delayed_without_wrap(t,x,tau):
    t=np.asarray(t,float);x=np.asarray(x,float)
    return np.interp(t-float(tau),t,x,left=0.0,right=0.0)

def estimate_series_windows(t,q,O,gas_series,heat,gas,window,compute_p=False,replicate=None,condition=None):
    rows=[];nin=window//3;maxlag=window//6
    for start in np.arange(0,t[-1]-window+1e-8,window):
        grid=start+2*np.arange(window//2);qx=np.interp(grid[:nin],t,q);gy=np.interp(grid,t,gas_series)
        if np.std(qx)<1e-7 or np.std(gy)<1e-7:continue
        x,y,phi=whiten(qx,gy)
        if np.std(x)<1e-9 or np.std(y)<1e-9:continue
        k,peak,signed=best_lag(x,y,maxlag,'pearson')
        if k is None:continue
        pval=np.nan;nnull=0
        if compute_p:
            null=[]
            for shift in range(len(y)):
                if shift<=maxlag or shift>=len(y)-maxlag:continue
                _,sc,_=best_lag(x,np.roll(y,-shift),maxlag,'pearson')
                if np.isfinite(sc):null.append(sc)
            nnull=len(null);pval=(1+sum(v>=peak for v in null))/(1+len(null)) if null else np.nan
        rows.append(dict(heat=heat,gas=gas,window_s=window,start_s=float(start),lag_s=float(2*k),
            signed_assoc=signed,peak_abs_assoc=peak,boundary_peak=bool(k==maxlag),zero_lag_peak=bool(k==0),
            p_value=pval,n_null=nnull,replicate=replicate,condition=condition))
    return rows

def prepare_simulation(raw,curves,outdir,seed=SEED,n_sim=N_SIM,rep_start=0):
    rng_master=np.random.default_rng(seed)
    raw=raw.sort_values(['attachment','row']);curve0=curves[curves.scenario=='H0']
    # Use the official H0 raw time grid and the frozen analytic oxygen input.
    base={}; trends={}; qflucs={}; residuals={}; sigma={}; block_seconds={}; block_samples={}; trend_diagnostics=[]
    for gas in GASES:
        residuals[gas]={}; efol=[]
        for heat,d in raw.groupby('heat',sort=False):
            d=d.sort_values('row'); t=(d.row.to_numpy(float)-1)*DT_H0
            gr=curve0[curve0.heat==heat].sort_values('time_s')
            q=gr.oxygen_flow_analytic_m3_s.to_numpy(float)
            y=d[gas].to_numpy(float);trend=rolling_trend(y,121);res=y-trend
            qtrend=rolling_trend(q,121);qf=q-qtrend
            residuals[gas][heat]=res-res.mean(); efol.append(acf_efold_seconds(res-res.mean(),DT_H0))
            base[(heat,gas)]={'t':t,'q':q,'O':gr.oxygen_cumulative_fit_m3.to_numpy(float),'raw_gas':y,'trend':trend,
                              'resid':res-res.mean(),'qtrend':qtrend,'qfluc':qf,'sigma':float(np.std(res-res.mean()))}
            trends[(heat,gas)]=trend;qflucs[(heat,gas)]=qf;sigma[(heat,gas)]=float(np.std(res-res.mean()))
        med=float(np.nanmedian(efol));secs=float(med);samples=max(2,int(np.ceil(secs/DT_H0)))
        block_seconds[gas]=secs;block_samples[gas]=samples
    for gas in GASES:
        for heat in HEATS:
            b=base[(heat,gas)];trend_diagnostics.append(dict(heat=heat,gas=gas,n=len(b['t']),raw_sd=np.std(b['raw_gas']),trend_sd=np.std(b['trend']),
                empirical_resid_sd=np.std(b['resid']),actual_residual_acf1=acf1(b['resid']),oxygen_fluct_sd=np.std(b['qfluc']),noise_block_acf_efold_s=block_seconds[gas],
                noise_block_samples=block_samples[gas],raw_missing_fraction=float(np.mean(~np.isfinite(b['raw_gas'])))))
    pd.DataFrame(trend_diagnostics).to_csv(outdir/'simulation_input_diagnostics.csv',index=False)

    # Deterministic noiseless recovery self-check at the raw observation layer.
    self_rows=[]
    for gas in GASES:
      for heat in HEATS:
        b=base[(heat,gas)];t=b['t'];qf=b['qfluc'];
        for tau in [0,10,20,30]:
          shifted=delayed_without_wrap(t,qf,tau)
          for direction,sgn in [('positive',1.),('negative',-1.)]:
            y=sgn*shifted; y=y+max(0.,-float(np.min(y)))+1e-4
            fit=k2_fit(t,y)(t)
            for window in [90,180]:
              wr=estimate_series_windows(t,b['q'],b['O'],fit,heat,gas,window,False)
              for z in wr:self_rows.append(dict(**z,true_lag_s=tau,direction=direction,scenario='noiseless',snr=np.inf))
    self_df=pd.DataFrame(self_rows);self_df.to_csv(outdir/'noiseless_self_check_windows.csv',index=False)
    if len(self_df):
        self_df['within_2s']=(np.abs(self_df.lag_s-self_df.true_lag_s)<=2).astype(float)
        ss=self_df.groupby(['gas','window_s','true_lag_s','direction'],as_index=False).agg(
            n_windows=('lag_s','size'),median_estimated_lag_s=('lag_s','median'),
            within_2s_rate=('within_2s','mean'),boundary_rate=('boundary_peak','mean'))
        cand=[]
        for window in [90,180]:
            ncan=sum(len(np.arange(0,base[(h,GASES[0])]['t'][-1]-window+1e-8,window)) for h in HEATS)
            cand.append((window,ncan))
        ss['candidate_windows']=ss.window_s.map(dict(cand));ss['estimated_window_rate']=ss.n_windows/ss.candidate_windows
        ss.to_csv(outdir/'noiseless_self_check_summary.csv',index=False)

    sim_rows=[];zero_rows=[];scale_rows=[]
    directions=[('positive',1.),('negative',-1.)]
    conditions=[(tau,direction,sgn,snr) for tau in [0,10,20,30] for direction,sgn in directions for snr in [1.,.5]]
    # Add a separate zero-response condition; response coefficient exactly zero.
    for gas in GASES:
      for rep in range(rep_start,rep_start+n_sim):
        rep_windows=[]
        # Zero-response control, all ten heats in the same replicate panel.
        for ih,heat in enumerate(HEATS):
            b=base[(heat,gas)];t=b['t'];qf=b['qfluc'];n=len(t)
            rng=np.random.default_rng(np.random.SeedSequence([seed,1 if gas=='CO' else 2,999,rep,ih]))
            noise=moving_block_noise(n,heat,residuals[gas],block_samples[gas],sigma[(heat,gas)],rng)
            yraw=b['trend']+noise; yraw=yraw+max(0.,-float(np.min(yraw)))+1e-4
            fit=k2_fit(t,yraw)(t)
            for window in [90,180]:
                wr=estimate_series_windows(t,b['q'],b['O'],fit,heat,gas,window,compute_p=True,replicate=rep,condition='zero_response')
                rep_windows.extend(wr)
        # BH per replicate x gas x window across all heat windows.
        for window in [90,180]:
            z=[r for r in rep_windows if r['window_s']==window]
            if z:
                qv=bh([r['p_value'] for r in z]);
                for rr,qv1 in zip(z,qv):rr['bh_q']=float(qv1);zero_rows.append(rr)
        # Injected response conditions.
        for ic,(tau,direction,sgn,snr) in enumerate(conditions):
          for ih,heat in enumerate(HEATS):
            b=base[(heat,gas)];t=b['t'];qf=b['qfluc'];n=len(t)
            rng=np.random.default_rng(np.random.SeedSequence([seed,1 if gas=='CO' else 2,ic+1,rep,ih]))
            noise=moving_block_noise(n,heat,residuals[gas],block_samples[gas],sigma[(heat,gas)],rng)
            shifted=delayed_without_wrap(t,qf,tau)
            valid=t>=tau
            qsd=float(np.std(shifted[valid])) if valid.any() else float(np.std(shifted))
            noise_sd=float(np.std(noise));beta=sgn*snr*noise_sd/max(qsd,1e-12)
            yraw=b['trend']+beta*shifted+noise
            yraw=yraw+max(0.,-float(np.min(yraw)))+1e-4
            fit=k2_fit(t,yraw)(t)
            for window in ([90,180] if tau==30 else [90]):
                wr=estimate_series_windows(t,b['q'],b['O'],fit,heat,gas,window,False,replicate=rep,
                    condition=f'tau{tau}_{direction}_snr{snr:g}')
                for z in wr:
                    z.update(true_lag_s=tau,direction=direction,snr=snr,response_beta=beta,signal_noise_ratio_target=snr,
                             signal_noise_ratio_realized=float(np.std(beta*shifted[valid])/max(noise_sd,1e-12)),
                             scenario='noisy_injected')
                    sim_rows.append(z)
            # Scale/ACF diagnostics for tau20 positive SNR=1 only.
            if tau==20 and direction=='positive' and snr==1.:
                syn_trend=rolling_trend(yraw,121);syn_resid=yraw-syn_trend
                scale_rows.append(dict(heat=heat,gas=gas,replicate=rep,synthetic_raw_sd=np.std(yraw),synthetic_trend_sd=np.std(syn_trend),
                    synthetic_residual_sd=np.std(syn_resid),synthetic_residual_acf1=acf1(syn_resid),synthetic_missing_fraction=float(np.mean(~np.isfinite(yraw))),
                    actual_raw_sd=np.std(b['raw_gas']),actual_trend_sd=np.std(b['trend']),actual_residual_sd=np.std(b['resid']),
                    actual_residual_acf1=acf1(b['resid']),actual_missing_fraction=float(np.mean(~np.isfinite(b['raw_gas']))),
                    injected_noise_sd=np.std(noise),injected_noise_acf1=acf1(noise),target_snr=snr,
                    realized_snr=float(np.std(beta*shifted[valid])/max(noise_sd,1e-12))))
        if rep%10==0: print(f'simulation {gas}: replicate {rep+1}/{rep_start+n_sim}',flush=True)
    sim=pd.DataFrame(sim_rows);zero=pd.DataFrame(zero_rows);scale=pd.DataFrame(scale_rows)
    sim.to_csv(outdir/'synthetic_recovery_windows.csv',index=False,float_format='%.8g')
    zero.to_csv(outdir/'zero_response_windows.csv',index=False,float_format='%.8g')
    scale.to_csv(outdir/'synthetic_scale_diagnostics_replicates.csv',index=False,float_format='%.8g')

    # Per-heat metrics first, then equal-weight across heats.
    if len(sim):
      sim['error_s']=sim.lag_s-sim.true_lag_s
      sim['direction_sign_correct']=np.where(sim.direction.eq('positive'),sim.signed_assoc>0,sim.signed_assoc<0).astype(float)
      ph=sim.groupby(['gas','window_s','true_lag_s','direction','snr','heat'],as_index=False).agg(
        n_windows=('lag_s','size'),bias_s=('error_s','mean'),mae_s=('error_s',lambda x:float(np.mean(np.abs(x)))),
        rmse_s=('error_s',lambda x:float(np.sqrt(np.mean(np.asarray(x)**2)))),within_2s_rate=('error_s',lambda x:float(np.mean(np.abs(x)<=2))),
        boundary_rate=('boundary_peak','mean'),upper_bound_rate=('boundary_peak','mean'),direction_sign_accuracy=('direction_sign_correct','mean'))
      # Candidate-window denominator is determined by the original heat time axis, not by successful estimates.
      cand=[]
      for heat in HEATS:
          tt=base[(heat, GASES[0])]['t']
          for window in [90,180]:
              ncan=len(np.arange(0,tt[-1]-window+1e-8,window))
              cand.append((heat,window,ncan))
      cand_df=pd.DataFrame(cand,columns=['heat','window_s','candidate_windows_per_replicate'])
      ph=ph.merge(cand_df,on=['heat','window_s'],how='left')
      ph['candidate_windows']=ph['candidate_windows_per_replicate']*n_sim
      ph['estimated_window_rate']=ph['n_windows']/ph['candidate_windows']
      ph.to_csv(outdir/'synthetic_recovery_by_heat.csv',index=False,float_format='%.8g')
      eq=ph.groupby(['gas','window_s','true_lag_s','direction','snr'],as_index=False).agg(
        n_heats=('heat','nunique'),equal_heat_bias_s=('bias_s','mean'),equal_heat_mae_s=('mae_s','mean'),
        equal_heat_rmse_s=('rmse_s','mean'),equal_heat_within_2s_rate=('within_2s_rate','mean'),
        equal_heat_boundary_rate=('boundary_rate','mean'),equal_heat_upper_bound_rate=('upper_bound_rate','mean'),
        mean_windows_per_heat=('n_windows','mean'))
      eq.to_csv(outdir/'synthetic_recovery_equal_heat.csv',index=False,float_format='%.8g')
    if len(zero):
      zh=zero.groupby(['gas','window_s','heat'],as_index=False).agg(n_windows=('lag_s','size'),median_peak_lag_s=('lag_s','median'),
          peak_q25_s=('lag_s',lambda x:float(x.quantile(.25))),peak_q75_s=('lag_s',lambda x:float(x.quantile(.75))),
          positive_peak_rate=('signed_assoc',lambda x:float(np.mean(np.asarray(x)>0))),negative_peak_rate=('signed_assoc',lambda x:float(np.mean(np.asarray(x)<0))),
          bh_false_positive_rate=('bh_q',lambda x:float(np.mean(np.asarray(x)<=.1))))
      zh.to_csv(outdir/'zero_response_by_heat.csv',index=False,float_format='%.8g')
      zeq=zh.groupby(['gas','window_s'],as_index=False).agg(n_heats=('heat','nunique'),equal_heat_median_peak_lag_s=('median_peak_lag_s','mean'),
          equal_heat_peak_q25_s=('peak_q25_s','mean'),equal_heat_peak_q75_s=('peak_q75_s','mean'),
          equal_heat_positive_peak_rate=('positive_peak_rate','mean'),equal_heat_negative_peak_rate=('negative_peak_rate','mean'),
          equal_heat_bh_false_positive_rate=('bh_false_positive_rate','mean'))
      zeq.to_csv(outdir/'zero_response_equal_heat.csv',index=False,float_format='%.8g')
      act=actual_summary(get_windows(curves,'H0',90,'pearson','old',True))
      a180=actual_summary(get_windows(curves,'H0',180,'pearson','old',True))
      act=pd.concat([act,a180],ignore_index=True)
      zact=[]
      for _,r in zeq.iterrows():
        ar=act[(act.gas==r.gas)&(act.window_s==r.window_s)]
        if len(ar):
          zact.append(dict(gas=r.gas,window_s=r.window_s,actual_median_peak_lag_s=float(ar.iloc[0].lag_median_s),
             zero_response_equal_heat_mean_peak_lag_s=float(r.equal_heat_median_peak_lag_s),
             zero_response_peak_q25_equal_heat_s=float(r.equal_heat_peak_q25_s),zero_response_peak_q75_equal_heat_s=float(r.equal_heat_peak_q75_s),
             zero_response_bh_false_positive_rate=float(r.equal_heat_bh_false_positive_rate),
             actual_positive_peak_rate=float(ar.iloc[0].positive_peak_n/ar.iloc[0].n_windows),
             actual_negative_peak_rate=float(ar.iloc[0].negative_peak_n/ar.iloc[0].n_windows)))
      pd.DataFrame(zact).to_csv(outdir/'actual_vs_zero_response.csv',index=False,float_format='%.8g')
    return sim,zero,self_df,scale,block_seconds,block_samples

def cdiff(t,y,half):
    a=np.maximum(t-half,t[0]);b=np.minimum(t+half,t[-1]);return (np.interp(b,t,y)-np.interp(a,t,y))/(b-a)

def round2_40(curves,outdir):
    methods_out=[]; actual_rows=[]
    for (scenario,heat),f in curves.groupby(['scenario','heat'],sort=False):
        f=f.sort_values('time_s');t=f.time_s.to_numpy(float);O=f.oxygen_cumulative_fit_m3.to_numpy(float)
        flow_methods={
          'analytic':f.oxygen_flow_analytic_m3_s.to_numpy(float),
          'raw_cum_30s':cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),15.),
          'raw_cum_60s':cdiff(t,f.oxygen_cumulative_raw_m3.to_numpy(float),30.),
          'fit_cum_30s':cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),15.),
          'fit_cum_60s':cdiff(t,f.oxygen_cumulative_fit_m3.to_numpy(float),30.)}
        for method,q in flow_methods.items():
          for win in WINDOWS:
            for offset in [0,win//2]:
              nin=win//3
              for gas in GASES:
                ygas=f[f'{gas}_flow_k2_m3_s'].to_numpy(float)
                for start in np.arange(offset,t[-1]-win+1e-8,win):
                    grid=start+2*np.arange(win//2);qr=np.interp(grid,t,q);yr=np.interp(grid,t,ygas)
                    if np.std(qr[:nin])<1e-7 or np.std(yr)<1e-7:continue
                    x,y,phi=whiten(qr[:nin],yr)
                    if np.std(x)<1e-9 or np.std(y)<1e-9:continue
                    k,score,signed=best_lag(x,y,win//6,'pearson')
                    if k is None:continue
                    progress=(np.interp(start+win/2,t,O)-O[0])/(O[-1]-O[0])
                    actual_rows.append(dict(kind='actual',scenario=scenario,heat=heat,method=method,gas=gas,window_s=win,offset_s=offset,
                        start_s=float(start),flow_m3_s=float(qr[:nin].mean()),oxygen_progress=float(progress),lag_s=float(2*k),
                        peak_abs_assoc=score,signed_assoc=signed,boundary_peak=bool(k==win//6)))
    w=pd.DataFrame(actual_rows);summ=[]
    for key,g in w.groupby(['scenario','method','gas','window_s','offset_s']):
        scenario,method,gas,win,offset=key;ds=[];parts=[]
        for heat,h in g.groupby('heat'):
            if len(h)<4:continue
            q1,q2=np.quantile(h.flow_m3_s,[1/3,2/3]);lo=h[h.flow_m3_s<=q1];hi=h[h.flow_m3_s>q2]
            if len(lo) and len(hi):ds.append((heat,float(hi.lag_s.median()-lo.lag_s.median())))
            x=h.flow_m3_s.to_numpy(float);y=h.lag_s.to_numpy(float);s=h.oxygen_progress.to_numpy(float)
            if len(x)<3:continue
            X=np.column_stack([np.ones(len(s)),s]);xr=x-X@np.linalg.lstsq(X,x,rcond=None)[0];yr=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
            parts.append((heat,float(xr@xr/len(xr)),float(xr@yr/len(yr))))
        d=np.array([z[1] for z in ds]);nz=d[d!=0];kneg=int((nz<0).sum());n=len(nz)
        signp=float(sum(math.comb(n,j) for j in range(kneg,n+1))/(2**n)) if n else np.nan
        den=sum(z[1] for z in parts);num=sum(z[2] for z in parts);slope=num/den if den>0 else np.nan
        loo=[]
        for omit,_,_ in parts:
            rest=[z for z in parts if z[0]!=omit];dd=sum(z[1] for z in rest)
            if dd>0:loo.append(sum(z[2] for z in rest)/dd)
        summ.append(dict(scenario=scenario,method=method,gas=gas,window_s=win,offset_s=offset,n_windows=len(g),n_heats=g.heat.nunique(),
            median_lag_s=g.lag_s.median(),lag_q25_s=g.lag_s.quantile(.25),lag_q75_s=g.lag_s.quantile(.75),
            median_high_minus_low_lag_s=float(np.median(d)) if len(d) else np.nan,n_heats_contrast=len(d),
            n_heats_high_flow_shorter=int((d<0).sum()) if len(d) else np.nan,n_nonzero_heat_contrasts=n,
            sign_test_one_sided_p=signp,stage_linear_adjusted_slope=slope,slope_loo_min=min(loo) if loo else np.nan,
            slope_loo_max=max(loo) if loo else np.nan,n_loo_negative=sum(z<0 for z in loo),n_loo=len(loo)))
    s=pd.DataFrame(summ)
    primary=s[s.offset_s==0].copy();pvals=primary.sign_test_one_sided_p.fillna(1.).to_numpy(float)
    order=np.argsort(pvals);adj=np.empty(len(pvals));run=0.
    for rank,idx in enumerate(order):run=max(run,(len(pvals)-rank)*pvals[idx]);adj[idx]=min(1.,run)
    primary['holm_p_40']=adj
    s=s.merge(primary[['scenario','method','gas','window_s','holm_p_40']],on=['scenario','method','gas','window_s'],how='left')
    s.loc[s.offset_s!=0,'holm_p_40']=np.nan
    w.to_csv(outdir/'round2_repaired_actual_windows.csv',index=False,float_format='%.8g')
    s.to_csv(outdir/'round2_repaired_actual_summary.csv',index=False,float_format='%.8g')
    return w,s

def make_figures(outdir,actual_windows,actual_summary,bins,match,sim_eq,zero_eq):
    figs=Path(outdir)/'figures';figs.mkdir(exist_ok=True)
    # Actual flow-progress overlap by heat, with low/high tercile tails.
    w=actual_windows[(actual_windows.scenario=='H0')&(actual_windows.window_s==90)&(actual_windows.metric=='pearson')].copy()
    fig,axs=plt.subplots(2,5,figsize=(17,7),sharex=True,sharey=True)
    for ax,heat in zip(axs.flat,HEATS):
        h=w[w.heat==heat];q1,q3=np.quantile(h.flow_m3_s,[.25,.75]);
        lab=np.where(h.flow_m3_s<q1,'low',np.where(h.flow_m3_s>q3,'high','middle'))
        for name,marker in [('low','o'),('middle','.'),('high','^')]:
            z=h[lab==name];ax.scatter(z.oxygen_progress_new,z.flow_m3_s,s=24,marker=marker,label=name,alpha=.8)
        ax.set_title(heat,fontsize=8);ax.grid(alpha=.2)
    axs[1,2].set_xlabel('Oxygen progress (new input-interval mean)');axs[0,0].set_ylabel('Mean oxygen flow (m³/s)');axs[1,0].set_ylabel('Mean oxygen flow (m³/s)')
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3);fig.suptitle('Q2.3 within-heat flow–progress overlap (H0, 90 s)',y=1.02)
    fig.tight_layout();fig.savefig(figs/'flow_progress_by_heat.png',dpi=180,bbox_inches='tight');plt.close(fig)
    # Actual lag vs zero-response null.
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for ax,gas in zip(axs,GASES):
        z=zero_eq[(zero_eq.gas==gas)&(zero_eq.window_s==90)]
        a=actual_summary[(actual_summary.scenario=='H0')&(actual_summary.window_s==90)&(actual_summary.metric=='pearson')&(actual_summary.gas==gas)]
        if len(z): ax.axvline(float(a.lag_median_s.iloc[0]),color='black',ls='--',label=f"actual median {a.lag_median_s.iloc[0]:.0f}s")
        ax.bar(['zero-response\nmedian peak','actual\nmedian peak'],[float(z.equal_heat_median_peak_lag_s.iloc[0]) if len(z) else np.nan,float(a.lag_median_s.iloc[0]) if len(a) else np.nan])
        ax.set_title(gas);ax.set_ylabel('Peak lag (s)');ax.grid(axis='y',alpha=.2);ax.legend()
    fig.suptitle('Actual peak lag compared with zero-response synthetic control');fig.tight_layout();fig.savefig(figs/'actual_vs_zero_response.png',dpi=180);plt.close(fig)
    # Recovery accuracy by SNR and window.
    if len(sim_eq):
        fig,axs=plt.subplots(1,2,figsize=(12,4),sharey=True)
        for ax,gas in zip(axs,GASES):
            z=sim_eq[(sim_eq.gas==gas)&(sim_eq.window_s==90)]
            for snr,gg in z.groupby('snr'):
                for direction,g2 in gg.groupby('direction'):
                    m=g2.groupby('true_lag_s').equal_heat_mae_s.mean()
                    ax.plot(m.index,m.values,marker='o',label=f'{direction}, SNR={snr:g}')
            ax.set_title(gas);ax.set_xlabel('Injected delay (s)');ax.grid(alpha=.2);ax.legend(fontsize=8)
        axs[0].set_ylabel('Equal-heat MAE (s)');fig.suptitle('Noisy semi-synthetic lag recovery, 90-s window');fig.tight_layout();fig.savefig(figs/'lag_recovery_mae.png',dpi=180);plt.close(fig)

def aggregate_simulation_batches(raw,curves,batch_dirs,outdir,n_sim,seed):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    files=['synthetic_recovery_windows.csv','zero_response_windows.csv','synthetic_scale_diagnostics_replicates.csv']
    for name in files:
        frames=[pd.read_csv(Path(b)/name) for b in batch_dirs if (Path(b)/name).exists()]
        pd.concat(frames,ignore_index=True).to_csv(outdir/name,index=False,float_format='%.8g')
    # Deterministic diagnostics and noiseless self-check are identical across batches; retain first.
    for name in ['simulation_input_diagnostics.csv','noiseless_self_check_windows.csv','noiseless_self_check_summary.csv']:
        src=Path(batch_dirs[0])/name
        if src.exists():
            import shutil;shutil.copyfile(src,outdir/name)
    sim=pd.read_csv(outdir/'synthetic_recovery_windows.csv')
    zero=pd.read_csv(outdir/'zero_response_windows.csv')
    scale=pd.read_csv(outdir/'synthetic_scale_diagnostics_replicates.csv')
    self_df=pd.read_csv(outdir/'noiseless_self_check_windows.csv')
    if len(self_df):
        self_df['within_2s']=(np.abs(self_df.lag_s-self_df.true_lag_s)<=2).astype(float)
        self_df['direction_sign_correct']=np.where(self_df.direction.eq('positive'),self_df.signed_assoc>0,self_df.signed_assoc<0).astype(float)
        ss=self_df.groupby(['gas','window_s','true_lag_s','direction'],as_index=False).agg(
            n_windows=('lag_s','size'),median_estimated_lag_s=('lag_s','median'),within_2s_rate=('within_2s','mean'),
            boundary_rate=('boundary_peak','mean'),direction_sign_accuracy=('direction_sign_correct','mean'))
        cand=[]
        for window in [90,180]:
            ncan=sum(len(np.arange(0,(raw[raw.heat==h].row.max()-1)*DT_H0-window+1e-8,window)) for h in HEATS)
            cand.append((window,ncan))
        ss['candidate_windows']=ss.window_s.map(dict(cand));ss['estimated_window_rate']=ss.n_windows/ss.candidate_windows
        ss.to_csv(outdir/'noiseless_self_check_summary.csv',index=False)
    scale.groupby(['heat','gas'],as_index=False).mean(numeric_only=True).to_csv(outdir/'synthetic_scale_diagnostics_summary.csv',index=False,float_format='%.8g')
    if len(sim):
        sim['error_s']=sim.lag_s-sim.true_lag_s
        sim['direction_sign_correct']=np.where(sim.direction.eq('positive'),sim.signed_assoc>0,sim.signed_assoc<0).astype(float)
        ph=sim.groupby(['gas','window_s','true_lag_s','direction','snr','heat'],as_index=False).agg(
            n_windows=('lag_s','size'),bias_s=('error_s','mean'),mae_s=('error_s',lambda x:float(np.mean(np.abs(x)))),
            rmse_s=('error_s',lambda x:float(np.sqrt(np.mean(np.asarray(x)**2)))),within_2s_rate=('error_s',lambda x:float(np.mean(np.abs(x)<=2))),
            boundary_rate=('boundary_peak','mean'),upper_bound_rate=('boundary_peak','mean'),direction_sign_accuracy=('direction_sign_correct','mean'))
        cand=[]
        for heat in HEATS:
            tt=(raw[raw.heat==heat].row.to_numpy(float)-1)*DT_H0
            for window in [90,180]:cand.append((heat,window,len(np.arange(0,tt[-1]-window+1e-8,window))))
        cand_df=pd.DataFrame(cand,columns=['heat','window_s','candidate_windows_per_replicate'])
        ph=ph.merge(cand_df,on=['heat','window_s'],how='left');ph['candidate_windows']=ph.candidate_windows_per_replicate*n_sim
        ph['estimated_window_rate']=ph.n_windows/ph.candidate_windows
        ph.to_csv(outdir/'synthetic_recovery_by_heat.csv',index=False,float_format='%.8g')
        eq=ph.groupby(['gas','window_s','true_lag_s','direction','snr'],as_index=False).agg(
            n_heats=('heat','nunique'),equal_heat_bias_s=('bias_s','mean'),equal_heat_mae_s=('mae_s','mean'),
            equal_heat_rmse_s=('rmse_s','mean'),equal_heat_within_2s_rate=('within_2s_rate','mean'),
            equal_heat_boundary_rate=('boundary_rate','mean'),equal_heat_upper_bound_rate=('upper_bound_rate','mean'),
            equal_heat_direction_sign_accuracy=('direction_sign_accuracy','mean'),equal_heat_estimated_window_rate=('estimated_window_rate','mean'),mean_windows_per_heat=('n_windows','mean'))
        eq.to_csv(outdir/'synthetic_recovery_equal_heat.csv',index=False,float_format='%.8g')
    if len(zero):
        zh=zero.groupby(['gas','window_s','heat'],as_index=False).agg(n_windows=('lag_s','size'),median_peak_lag_s=('lag_s','median'),
            peak_q25_s=('lag_s',lambda x:float(x.quantile(.25))),peak_q75_s=('lag_s',lambda x:float(x.quantile(.75))),
            positive_peak_rate=('signed_assoc',lambda x:float(np.mean(np.asarray(x)>0))),negative_peak_rate=('signed_assoc',lambda x:float(np.mean(np.asarray(x)<0))),
            bh_false_positive_rate=('bh_q',lambda x:float(np.mean(np.asarray(x)<=.1))))
        zh.to_csv(outdir/'zero_response_by_heat.csv',index=False,float_format='%.8g')
        zeq=zh.groupby(['gas','window_s'],as_index=False).agg(n_heats=('heat','nunique'),equal_heat_median_peak_lag_s=('median_peak_lag_s','mean'),
            equal_heat_peak_q25_s=('peak_q25_s','mean'),equal_heat_peak_q75_s=('peak_q75_s','mean'),
            equal_heat_positive_peak_rate=('positive_peak_rate','mean'),equal_heat_negative_peak_rate=('negative_peak_rate','mean'),
            equal_heat_bh_false_positive_rate=('bh_false_positive_rate','mean'))
        zeq.to_csv(outdir/'zero_response_equal_heat.csv',index=False,float_format='%.8g')
        act=pd.concat([actual_summary(get_windows(curves,'H0',90,'pearson','old',True)),actual_summary(get_windows(curves,'H0',180,'pearson','old',True))],ignore_index=True)
        rows=[]
        for _,r in zeq.iterrows():
            ar=act[(act.gas==r.gas)&(act.window_s==r.window_s)]
            if len(ar):rows.append(dict(gas=r.gas,window_s=r.window_s,actual_median_peak_lag_s=float(ar.iloc[0].lag_median_s),
                zero_response_equal_heat_mean_peak_lag_s=float(r.equal_heat_median_peak_lag_s),
                zero_response_peak_q25_equal_heat_s=float(r.equal_heat_peak_q25_s),zero_response_peak_q75_equal_heat_s=float(r.equal_heat_peak_q75_s),
                zero_response_bh_false_positive_rate=float(r.equal_heat_bh_false_positive_rate),
                actual_positive_peak_rate=float(ar.iloc[0].positive_peak_n/ar.iloc[0].n_windows),
                actual_negative_peak_rate=float(ar.iloc[0].negative_peak_n/ar.iloc[0].n_windows)))
        pd.DataFrame(rows).to_csv(outdir/'actual_vs_zero_response.csv',index=False,float_format='%.8g')
        zero['reject_q_le_0_1']=zero.bh_q<=.1
        fpr=[]
        for (gas,window),g in zero.groupby(['gas','window_s']):
            rep=g.groupby('replicate').reject_q_le_0_1.agg(['sum','any','count'])
            fpr.append(dict(gas=gas,window_s=window,n_replicates=len(rep),replicates_with_any_discovery=int(rep['any'].sum()),
                replicates_with_any_q_rejection_rate=float(rep['any'].mean()),mean_q_rejections_per_replicate=float(rep['sum'].mean()),
                mean_window_q_rejection_rate=float((rep['sum']/rep['count']).mean()),
                max_q_rejections_in_one_replicate=int(rep['sum'].max())))
        pd.DataFrame(fpr).to_csv(outdir/'zero_response_false_positive_summary.csv',index=False,float_format='%.8g')
        cond=sim.groupby(['gas','true_lag_s','direction','snr','window_s']).replicate.nunique().reset_index(name='replicates')
        cond.to_csv(outdir/'simulation_condition_coverage.csv',index=False)
    return sim,zero,scale

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--process',type=Path,required=True);ap.add_argument('--curves',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--archive',type=Path,default=Path('/mnt/data/Data.zip'));ap.add_argument('--simulations',type=int,default=N_SIM);ap.add_argument('--seed',type=int,default=SEED)
    a=ap.parse_args();out=a.out;out.mkdir(parents=True,exist_ok=True)
    raw=pd.read_csv(a.process);curves=pd.read_csv(a.curves)
    # Input and key uniqueness checks.
    process_key=raw.duplicated(['heat','row']).sum();curve_key=curves.duplicated(['scenario','heat','row']).sum()
    if len(raw)!=13589 or raw.heat.nunique()!=10 or process_key:raise ValueError('Raw process input failed row/heat/key checks')
    if len(curves)!=27178 or curves.heat.nunique()!=10 or curves.scenario.nunique()!=2 or curve_key:raise ValueError('Reconstructed curves failed row/heat/scenario/key checks')
    cfg=CONFIG.copy();cfg['simulation_replicates']=a.simulations;cfg['random_seed']=a.seed
    (out/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    # Actual windows: all official timebase/window/metric combinations; progress old/new side-by-side.
    actual_parts=[]
    actual_specs=[('H0',90,'pearson',True),('H0',180,'pearson',True),('H1',90,'pearson',True),
                  ('H1',180,'pearson',True),('H0',90,'spearman',False)]
    for scen,win,metric,do_p in actual_specs:
        wo=get_windows(curves,scen,win,metric,'old',do_p);wn=get_windows(curves,scen,win,metric,'new',False)
        if len(wo):
            keys=['metric','scenario','window_s','heat','gas','start_s']
            merged=wo.merge(wn[keys+['oxygen_progress']],on=keys,suffixes=('_old','_new'))
            merged=merged.rename(columns={'oxygen_progress_old':'oxygen_progress_old','oxygen_progress_new':'oxygen_progress_new'})
            actual_parts.append(merged)
    actual=pd.concat(actual_parts,ignore_index=True);print('ACTUAL WINDOWS DONE',len(actual),flush=True);actual.to_csv(out/'actual_lag_windows.csv',index=False,float_format='%.10g')
    summary=actual_summary(actual.rename(columns={'oxygen_progress_old':'oxygen_progress'}));summary.to_csv(out/'actual_lag_summary.csv',index=False,float_format='%.10g')
    bins=flow_bins_table(actual.rename(columns={'oxygen_progress_old':'oxygen_progress'}));print('ACTUAL SUMMARIES DONE',len(summary),len(bins),flush=True);bins.to_csv(out/'actual_flow_bins.csv',index=False,float_format='%.10g')
    # Flow matching and overlap; old/new progress, calipers 0.10 and 0.15.
    match_rows=[];match_summ=[];overlap_rows=[]
    for progress_def,col in [('old','oxygen_progress_old'),('new','oxygen_progress_new')]:
      g=actual[(actual.scenario=='H0')&(actual.window_s==90)&(actual.metric=='pearson')].copy();g['oxygen_progress']=g[col]
      overlap=overlap_diagnostics(g);overlap['progress_definition']=progress_def;overlap_rows.append(overlap)
      for gas in GASES:
        gg=g[g.gas==gas]
        for cal in [.10,.15]:
          pp=match_pairs(gg,cal,'oxygen_progress')
          if len(pp):
            pp['progress_definition']=progress_def;match_rows.append(pp)
            heat_eff=pp.groupby('heat').lag_diff_high_minus_low_s.median().to_numpy(float)
            rng=np.random.default_rng(a.seed+int(cal*1000)+(0 if gas=='CO' else 100)+(0 if progress_def=='old' else 1000))
            boot=np.median(heat_eff[rng.integers(0,len(heat_eff),size=(20000,len(heat_eff)))],axis=1)
            nz=heat_eff[np.abs(heat_eff)>1e-12];neg=int(np.sum(nz<0));n=len(nz)
            sp=binomtest(neg,n,.5,alternative='greater').pvalue if n else np.nan
            match_summ.append(dict(progress_definition=progress_def,gas=gas,caliper=cal,n_pairs=len(pp),n_heats=pp.heat.nunique(),
                median_progress_gap=pp.progress_gap.median(),max_progress_gap=pp.progress_gap.max(),median_flow_diff_m3_s=pp.flow_diff_m3_s.median(),
                median_pair_lag_diff_s=pp.lag_diff_high_minus_low_s.median(),median_heat_lag_diff_s=np.median(heat_eff),
                bootstrap_heat_median_ci_low=np.quantile(boot,.025),bootstrap_heat_median_ci_high=np.quantile(boot,.975),
                n_heats_high_flow_shorter=int((heat_eff<0).sum()),n_heats_high_flow_longer=int((heat_eff>0).sum()),
                sign_test_p_one_sided=sp))
          else:
            match_summ.append(dict(progress_definition=progress_def,gas=gas,caliper=cal,n_pairs=0,n_heats=0))
    if match_rows:pd.concat(match_rows,ignore_index=True).to_csv(out/'matched_pairs_old_new.csv',index=False,float_format='%.10g')
    pd.DataFrame(match_summ).to_csv(out/'matching_summary_old_new.csv',index=False,float_format='%.10g')
    print('MATCHING DONE',len(match_summ),flush=True)
    ov=pd.concat(overlap_rows,ignore_index=True)
    ov.to_csv(out/'flow_progress_overlap_by_heat.csv',index=False,float_format='%.10g')
    ovsum=ov.groupby('progress_definition',as_index=False).agg(n_heats=('heat','nunique'),heats_with_progress_overlap=('any_progress_overlap','sum'),
        total_low_windows_inside_overlap=('low_windows_inside_overlap','sum'),total_high_windows_inside_overlap=('high_windows_inside_overlap','sum'),
        median_overlap_progress_width=('overlap_progress_width','median'))
    ovsum.to_csv(out/'flow_progress_overlap_summary.csv',index=False,float_format='%.10g')
    # Fourth-round original regression specifications, old and corrected progress.
    reg_rows=[];heat_slope_rows=[]
    for progress_def,col in [('old','oxygen_progress_old'),('new','oxygen_progress_new')]:
      for gas in GASES:
        g=actual[(actual.scenario=='H0')&(actual.window_s==90)&(actual.metric=='pearson')&(actual.gas==gas)].copy()
        g['oxygen_progress']=g[col]
        for degree in [1,2]:
          for boundary in ['all','interior']:
            r=pooled_cluster_reg(g,degree,boundary)
            if r:reg_rows.append(dict(progress_definition=progress_def,gas=gas,**r))
          hs=heat_slopes(g,degree)
          if len(hs):
            hs['progress_definition']=progress_def;hs['gas']=gas;hs['degree']=degree;heat_slope_rows.append(hs)
    reg=pd.DataFrame(reg_rows);print('REGRESSION DONE',len(reg),flush=True);reg.to_csv(out/'flow_lag_regression_old_new.csv',index=False,float_format='%.10g')
    if heat_slope_rows:pd.concat(heat_slope_rows,ignore_index=True).to_csv(out/'heat_level_flow_lag_slopes_old_new.csv',index=False,float_format='%.10g')
    # Second-round 40 primary tests with the syntax-only fix; save complete family and Holm correction.
    _,r2=round2_40(curves,out);print('ROUND2 DONE',len(r2),flush=True);r2.to_csv(out/'round2_repaired_40_test_family.csv',index=False,float_format='%.10g')
    # Simulation and controls.
    print('SIM START',flush=True)
    batch_dirs=[];batch_logs=[]
    if a.simulations>=1:
        remaining=a.simulations;start=0;bi=0
        while remaining>0:
            take=min(10,remaining);bd=out/'simulation_batches'/f'batch_{bi+1:02d}';bd.mkdir(parents=True,exist_ok=True)
            prepare_simulation(raw,curves,bd,a.seed,take,start);batch_dirs.append(bd)
            batch_logs.append(dict(batch=bi+1,replicate_start=start,replicate_end=start+take-1,n_replicates=take,seed=a.seed))
            remaining-=take;start+=take;bi+=1
        sim,zero,scale=aggregate_simulation_batches(raw,curves,batch_dirs,out,a.simulations,a.seed)
        pd.DataFrame(batch_logs).to_csv(out/'simulation_batch_log.csv',index=False)
        block_seconds={g:float(pd.read_csv(out/'simulation_input_diagnostics.csv').query('gas==@g').noise_block_acf_efold_s.median()) for g in GASES}
        block_samples={g:int(pd.read_csv(out/'simulation_input_diagnostics.csv').query('gas==@g').noise_block_samples.median()) for g in GASES}
        self_df=pd.read_csv(out/'noiseless_self_check_windows.csv')
        import shutil;shutil.rmtree(out/'simulation_batches',ignore_errors=True)
    else:
        raise ValueError('simulation count must be >=1 for the locked final validation')
    print('SIM DONE',len(sim),len(zero),flush=True)
    cfg=json.loads((out/'config.json').read_text())
    diag=pd.read_csv(out/'simulation_input_diagnostics.csv')
    cfg['noise_block_seconds']={g:float(diag.loc[diag.gas==g,'noise_block_acf_efold_s'].median()) for g in GASES}
    cfg['noise_block_samples']={g:int(diag.loc[diag.gas==g,'noise_block_samples'].median()) for g in GASES}
    cfg['simulation_replicates']=a.simulations;cfg['random_seed']=a.seed
    cfg['response_amplitude']='beta = sign * target_SNR * sd(noise) / sd(shifted_oxygen_fluctuation over t>=tau); realized raw-layer SNR equals target by construction'
    (out/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    # Environment and checks.
    (out/'environment.txt').write_text(f'Python {platform.python_version()}\nnumpy {np.__version__}\npandas {pd.__version__}\nscipy {scipy.__version__}\nmatplotlib {matplotlib.__version__}\n',encoding='utf-8')
    process_bytes=Path(a.process).read_bytes()
    git_blob_sha1=hashlib.sha1(b'blob '+str(len(process_bytes)).encode()+b'\0'+process_bytes).hexdigest()
    checks={
      'raw_rows':len(raw),'raw_heats':int(raw.heat.nunique()),'raw_unique_heat_row_key':bool(process_key==0),
      'process_csv_github_blob_sha1':git_blob_sha1,'process_csv_matches_github_main_blob':bool(git_blob_sha1=='25a8f164274ff96e868ff447d28a4d1c80353b50'),
      'curve_rows':len(curves),'curve_heats':int(curves.heat.nunique()),'curve_scenarios':sorted(curves.scenario.unique()),
      'curve_unique_scenario_heat_row_key':bool(curve_key==0),'raw_missing_values':int(raw.isna().sum().sum()),
      'curve_missing_values':int(curves.isna().sum().sum()),'curve_sha256':sha256(a.curves),'process_sha256':sha256(a.process),
      'source_data_zip_sha256':sha256(a.archive),'noise_block_seconds':block_seconds,'noise_block_samples':block_samples,
      'round2_primary_tests':int(((r2.offset_s==0)).sum()),'round2_primary_holm_tests':int(r2.loc[r2.offset_s==0,'holm_p_40'].notna().sum()),
      'simulation_replicates_per_condition':a.simulations,'seed':a.seed}
    (out/'qa_manifest.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    # Config / diagnostics summary.
    sim_eq=pd.read_csv(out/'synthetic_recovery_equal_heat.csv') if (out/'synthetic_recovery_equal_heat.csv').exists() else pd.DataFrame()
    zero_eq=pd.read_csv(out/'zero_response_equal_heat.csv') if (out/'zero_response_equal_heat.csv').exists() else pd.DataFrame()
    print('FIGURES START',flush=True);make_figures(out,actual,summary,bins,pd.DataFrame(match_summ),sim_eq,zero_eq);print('FIGURES DONE',flush=True)
    # Difference audit for hashes / archived input references.
    diff=pd.DataFrame([
      dict(item='Data.zip SHA256',expected='cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727',observed=checks['source_data_zip_sha256'],status='MATCH' if checks['source_data_zip_sha256']=='cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727' else 'MISMATCH'),
      dict(item='Q2.3 Round3 input curves SHA256',expected='c89a8e34dc5f83385d6c06d356d2821f53e67b63197aa00b6dca9cf0f837f306',observed=checks['curve_sha256'],status='MATCH' if checks['curve_sha256']=='c89a8e34dc5f83385d6c06d356d2821f53e67b63197aa00b6dca9cf0f837f306' else 'MISMATCH'),
      dict(item='Q2 relationship study input curves SHA256',expected='9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49',observed=checks['curve_sha256'],status='HISTORICAL_INPUT_NOT_AVAILABLE'),
      dict(item='Q2 Round4 input curves SHA256',expected='f7280a7ee7276c6281b2c5f85b9fb003e73ea51ad8473656873a95f07fd016b4',observed=checks['curve_sha256'],status='HISTORICAL_INPUT_NOT_AVAILABLE'),
      dict(item='Q2.1 manifest process_csv_sha256',expected='f82e726b03f8eceb148a25440f34fe8385af28e95a7f551149b513ecf560c789',observed=checks['process_sha256'],status='MISMATCH_MANIFEST_STALE_OR_OTHER_EXPORT'),
      dict(item='GitHub main process.csv blob SHA1',expected='25a8f164274ff96e868ff447d28a4d1c80353b50',observed=checks['process_csv_github_blob_sha1'],status='MATCH' if checks['process_csv_matches_github_main_blob'] else 'MISMATCH'),
    ])
    diff.to_csv(out/'historical_input_hash_audit.csv',index=False)
    print('FINAL VALIDATION CORE COMPLETE')
    print('Actual lag summary H0/90/Pearson:')
    print(summary[(summary.scenario=='H0')&(summary.window_s==90)&(summary.metric=='pearson')].to_string(index=False))
    print('Matching summary:');print(pd.DataFrame(match_summ).to_string(index=False))
    print('Regression:');print(reg.to_string(index=False))
    print('Simulation equal heat:');print(sim_eq[(sim_eq.window_s==90)].to_string(index=False) if len(sim_eq) else 'none')
    print('Zero response equal heat:');print(zero_eq.to_string(index=False) if len(zero_eq) else 'none')

if __name__=='__main__':main()
