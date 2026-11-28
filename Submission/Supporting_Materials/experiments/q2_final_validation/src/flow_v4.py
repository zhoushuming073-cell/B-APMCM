from __future__ import annotations
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import lsq_linear
from scipy.integrate import trapezoid

def fit_integral_bspline(O,t,spacing_s,penalty,train_mask,irls_iter=2):
    O=np.asarray(O,dtype=float); t=np.asarray(t,dtype=float); train_mask=np.asarray(train_mask,dtype=bool)
    knots=np.r_[[t[0]]*4,np.arange(t[0]+spacing_s,t[-1],spacing_s),[t[-1]]*4]
    n_basis=len(knots)-4; basis=BSpline(knots,np.eye(n_basis),3,extrapolate=False)
    B=np.asarray(basis(t)); anti=basis.antiderivative(); A=np.asarray(anti(t)-anti(t[0])); D2=np.diff(np.eye(n_basis),n=2,axis=0)
    weights=np.ones(train_mask.sum()); coef=np.zeros(n_basis); success=False
    for _ in range(irls_iter):
        X=np.vstack([A[train_mask]*np.sqrt(weights)[:,None],np.sqrt(penalty)*D2])
        y=np.r_[(O[train_mask]-O[0])*np.sqrt(weights),np.zeros(len(D2))]
        fit=lsq_linear(X,y,bounds=(0,np.inf),max_iter=120); coef=fit.x; success=bool(fit.success)
        residual=O[train_mask]-O[0]-A[train_mask]@coef; center=np.median(residual)
        scale=max(1.4826*np.median(np.abs(residual-center)),0.1)
        weights=np.minimum(1.,1.345*scale/np.maximum(np.abs(residual-center),1e-9))
    q=B@coef; O_hat=O[0]+A@coef
    return q,O_hat,success

def make_blocks(n):
    blocks=[]
    for frac in (0.12,0.27,0.42,0.57,0.72,0.87):
        start=max(20,int(n*frac)); end=min(n-20,start+16); blocks.append(np.arange(start,end))
    return blocks

def cv_grid(O,t,base_mask,folds,grid):
    rows=[]
    for spacing_s,penalty in grid:
        fold_mae=[]
        for held in folds:
            train_mask=base_mask.copy(); train_mask[held]=False
            _,O_hat,_=fit_integral_bspline(O,t,spacing_s,penalty,train_mask,irls_iter=2)
            fold_mae.append(float(np.mean(np.abs(O_hat[held]-O[held]))))
        rows.append(dict(spacing_s=spacing_s,penalty=penalty,cv_MAE_m3=float(np.mean(fold_mae)),cv_MAE_se_m3=float(np.std(fold_mae,ddof=1)/np.sqrt(len(fold_mae))),**{f'fold_{i+1}_MAE_m3':v for i,v in enumerate(fold_mae)}))
    return pd.DataFrame(rows).sort_values(['cv_MAE_m3','spacing_s','penalty']).reset_index(drop=True)

def run(input_path:Path,output_dir:Path):
    output_dir.mkdir(parents=True,exist_ok=True); frame=pd.read_csv(input_path).sort_values(['attachment','row'])
    required={'row','O','CO2','CO','attachment','heat'}; missing=required-set(frame.columns)
    if missing: raise ValueError(f'Missing required columns: {sorted(missing)}')
    heats=list(frame.heat.drop_duplicates())
    if len(heats)!=10 or heats[4]!='6C02057': raise ValueError(f'Unexpected heat order/count: {heats}')
    grid=[(s,lam) for s in (5.,10.,20.) for lam in (.1,1.,10.,100.)]
    profile_rows=[]; summary_rows=[]; sensitivity_rows=[]; cv_rows=[]; validation_rows=[]
    for scenario in ('H0','H1'):
      for heat in heats:
        d=frame[frame.heat==heat].sort_values('row'); O=d.O.to_numpy(float); row=d.row.to_numpy(float)
        dt=1. if scenario=='H1' and heat=='6C02057' else .5; t=(row-1.)*dt; n=len(O); blocks=make_blocks(n)
        full_mask=np.ones(n,dtype=bool); folds=[blocks[k] for k in (0,2,3,5)]; full_cv=cv_grid(O,t,full_mask,folds,grid)
        best=full_cv.iloc[0]; cutoff=float(best.cv_MAE_m3+best.cv_MAE_se_m3); near=full_cv[full_cv.cv_MAE_m3<=cutoff].copy()
        if near.empty: near=full_cv.iloc[[0]].copy()
        spacing_best=float(best.spacing_s); penalty_best=float(best.penalty)
        q_best,O_best,success=fit_integral_bspline(O,t,spacing_best,penalty_best,full_mask,irls_iter=2)
        q_models=[]
        for _,candidate in near.iterrows():
            qc,_,_=fit_integral_bspline(O,t,float(candidate.spacing_s),float(candidate.penalty),full_mask,irls_iter=2); q_models.append(qc)
            sensitivity_rows.append(dict(scenario=scenario,heat=heat,spacing_s=float(candidate.spacing_s),penalty=float(candidate.penalty),cv_MAE_m3=float(candidate.cv_MAE_m3),cv_MAE_se_m3=float(candidate.cv_MAE_se_m3),within_one_SE=True,is_best=float(candidate.spacing_s)==spacing_best and float(candidate.penalty)==penalty_best))
        q_models=np.asarray(q_models); q_median=np.median(q_models,axis=0); q_p10=np.quantile(q_models,.1,axis=0); q_p90=np.quantile(q_models,.9,axis=0)
        def secant_rate(seconds):
            lag=max(1,int(round(seconds/dt))); left=np.maximum(0,np.arange(n)-lag); return (O-O[left])/np.maximum(t-t[left],dt)
        q10,q30,q60=secant_rate(10),secant_rate(30),secant_rate(60)
        for i in range(n):
            profile_rows.append(dict(scenario=scenario,heat=heat,row=int(row[i]),time_s=float(t[i]),cumulative_O_raw_m3=float(O[i]),cumulative_O_fit_m3=float(O_best[i]),flow_analytic_best_m3_s=float(q_best[i]),flow_sensitivity_median_m3_s=float(q_median[i]),flow_sensitivity_p10_m3_s=float(q_p10[i]),flow_sensitivity_p90_m3_s=float(q_p90[i]),observed_secant_10s_m3_s=float(q10[i]),observed_secant_30s_m3_s=float(q30[i]),observed_secant_60s_m3_s=float(q60[i])))
        integral_q=float(trapezoid(q_best,t))
        summary_rows.append(dict(scenario=scenario,heat=heat,dt_s_per_row=dt,n_rows=n,selected_spacing_s=spacing_best,selected_penalty=penalty_best,full_CV_MAE_m3=float(best.cv_MAE_m3),full_CV_MAE_SE_m3=float(best.cv_MAE_se_m3),near_optimal_model_count=int(len(near)),flow_min_m3_s=float(np.min(q_best)),flow_median_m3_s=float(np.median(q_best)),flow_mean_m3_s=float(np.mean(q_best)),flow_p95_m3_s=float(np.quantile(q_best,.95)),flow_max_m3_s=float(np.max(q_best)),integral_q_trapezoid_m3=integral_q,fitted_cumulative_delta_m3=float(O_best[-1]-O_best[0]),numeric_mass_balance_error_m3=integral_q-float(O_best[-1]-O_best[0]),mean_one_SE_sensitivity_band_width_m3_s=float(np.mean(q_p90-q_p10)),fit_success=success))
        me=abs(integral_q-float(O_best[-1]-O_best[0])); mono=int(np.sum(np.diff(O_best)<-1e-10)); neg=int(np.sum(q_best<-1e-10))
        validation_rows.append(dict(scenario=scenario,heat=heat,min_analytic_flow=float(np.min(q_best)),negative_flow_points=neg,cumulative_monotonic_violations=mono,abs_numeric_mass_balance_error_m3=me,fit_success=bool(success),pass_nonnegative=bool(np.min(q_best)>=-1e-10),pass_monotonic=bool(mono==0),pass_mass_balance_0_1m3=bool(me<=.1)))
        for _,r in full_cv.iterrows(): cv_rows.append(dict(scenario=scenario,heat=heat,spacing_s=float(r.spacing_s),penalty=float(r.penalty),cv_MAE_m3=float(r.cv_MAE_m3),cv_MAE_se_m3=float(r.cv_MAE_se_m3)))
    pd.DataFrame(profile_rows).to_csv(output_dir/'flow_profile.csv',index=False)
    pd.DataFrame(summary_rows).to_csv(output_dir/'flow_summary.csv',index=False)
    pd.DataFrame(sensitivity_rows).to_csv(output_dir/'flow_one_se_sensitivity_set.csv',index=False)
    pd.DataFrame(cv_rows).to_csv(output_dir/'flow_cv_grid.csv',index=False)
    pd.DataFrame(validation_rows).to_csv(output_dir/'flow_validation_checks.csv',index=False)
    return pd.DataFrame(profile_rows),pd.DataFrame(summary_rows),pd.DataFrame(validation_rows)
