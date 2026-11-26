#!/usr/bin/env python3
"""Matched-row LOHO finite distributed-lag models and flow-conditioned delay models."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
from scipy.optimize import lsq_linear
from scipy.stats import spearmanr,pearsonr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

p=argparse.ArgumentParser(); p.add_argument('--cache-dir',default='cache'); p.add_argument('--out',default='results'); a=p.parse_args()
cache=Path(a.cache_dir).resolve(); out=Path(a.out).resolve(); tables=out/'tables'; tables.mkdir(parents=True,exist_ok=True)
D={}
for f in sorted(cache.glob('heats_*.npz')):
 z=np.load(f); h=f.stem.replace('heats_',''); D[h]={'seq':z['seq'],'proc':z['proc']}
if len(D)!=10: raise ValueError(f'Expected 10 caches, found {len(D)}')
TARGET='6C02057'; lagdf=pd.read_csv(tables/'lag_estimates.csv')
# Cross-heat fixed, unconstrained flow-conditioned, and non-increasing flow-conditioned lag models.
lag_models=[]; coeffs=[]
for gas in ['CO2','CO']:
 for scenario,clock in [('H0','H0'),('H1','H1')]:
  other=lagdf[(lagdf.heat!=TARGET)&(lagdf.heat!='6B03205')&(lagdf.clock=='official')&(lagdf.gas==gas)]
  fifth=lagdf[(lagdf.heat==TARGET)&(lagdf.clock==clock)&(lagdf.gas==gas)]
  s=pd.concat([other,fifth],ignore_index=True).dropna(subset=['lag_s','flow_median_m3_s'])
  x=s.flow_median_m3_s.to_numpy(); y=s.lag_s.to_numpy(); X=np.c_[np.ones(len(x)),x]
  b=np.linalg.lstsq(X,y,rcond=None)[0]; bm=lsq_linear(X,y,bounds=([-np.inf,-np.inf],[np.inf,0.])).x
  sr=spearmanr(x,y); pr=pearsonr(x,y)
  coeffs.append({'scenario':scenario,'gas':gas,'n_heats':len(y),'fixed_mean_lag_s':np.mean(y),'fixed_median_lag_s':np.median(y),
   'target_lag_s':float(fifth.lag_s.iloc[0]),'target_corr':float(fifth['corr'].iloc[0]),'target_shift_p':float(fifth.shift_p.iloc[0]),
   'linear_intercept':b[0],'linear_slope_s_per_m3s':b[1],'monotone_intercept':bm[0],'monotone_slope_s_per_m3s':bm[1],
   'spearman_rho':sr.statistic,'spearman_p':sr.pvalue,'pearson_r':pr.statistic,'pearson_p':pr.pvalue})
  for model in ['fixed_mean','linear_unconstrained','linear_nonincreasing']:
   errs=[]
   for i in range(len(y)):
    m=np.arange(len(y))!=i; xx=x[m]; yy=y[m]
    if model=='fixed_mean': pred=float(np.mean(yy))
    elif model=='linear_unconstrained':
     bb=np.linalg.lstsq(np.c_[np.ones(len(xx)),xx],yy,rcond=None)[0]; pred=float(bb[0]+bb[1]*x[i])
    else:
     bb=lsq_linear(np.c_[np.ones(len(xx)),xx],yy,bounds=([-np.inf,-np.inf],[np.inf,0.])).x; pred=float(bb[0]+bb[1]*x[i])
    errs.append(abs(pred-y[i]))
   lag_models.append({'scenario':scenario,'gas':gas,'model':model,'n_heats':len(y),'LOHO_MAE_s':float(np.mean(errs)),
    'spearman_rho':float(sr.statistic),'spearman_p':float(sr.pvalue),'pearson_r':float(pr.statistic),'pearson_p':float(pr.pvalue)})
pd.DataFrame(coeffs).to_csv(tables/'lag_model_coefficients_and_summary.csv',index=False)
pd.DataFrame(lag_models).to_csv(tables/'flow_conditioned_lag_models.csv',index=False)

# Model features use native rows. H0/H1 keep identical raw-row train/test masks; only the physical lag-to-row mapping changes.
def lagged(x,k):
 n=len(x); z=np.full(n,np.nan); k=int(k)
 if k==0:return np.asarray(x,float).copy()
 if k<n:z[k:]=x[:-k]
 return z

def build(scenario):
 out={}
 for h,d in D.items():
  dt=1. if h==TARGET and scenario=='H1' else .5
  P=d['proc']; O=P[:,0]; G=P[:,1:3]; n=len(O)
  # causal 10-second oxygen-flow estimate
  from analysis_helpers import flow_slope
  q=flow_slope(O,dt,10.,True)
  lr=lambda s:int(round(s/dt))
  persistence=np.column_stack([lagged(G[:,0],lr(s)) for s in [1.,5.,10.]]+
                              [lagged(G[:,1],lr(s)) for s in [1.,5.,10.]])
  qlag=np.column_stack([lagged(q,lr(s)) for s in [0.,1.,2.,5.,10.,20.,30.,40.]])
  hist=[]
  for aa,bb in [(0,5),(5,10),(10,20),(20,40)]:
   k1,k2=lr(aa),lr(bb); v=np.full(n,np.nan)
   for i in range(n):
    lo=max(0,i-k2); hi=max(0,i-k1)
    if hi>lo:v[i]=np.mean(q[lo:hi])
   hist.append(v)
  hist=np.column_stack(hist); valid=np.arange(n)>=83
  features={'M0_persistence':persistence,'M1_current_flow':np.column_stack([persistence,qlag[:,0]]),
   'M2_distributed_lag':np.column_stack([persistence,qlag]),'M3_history_windows':np.column_stack([persistence,hist])}
  out[h]={'heat':h,'features':{k:v[valid] for k,v in features.items()},'y':G[valid]}
 return out

def fit_eval(train,test,feature,gi,alphas=(.1,10.)):
 scores=[]
 for alpha in alphas:
  vals=[]
  for vh in train:
   tr=[p for p in train if p['heat']!=vh['heat']]
   X=np.vstack([p['features'][feature] for p in tr]); y=np.concatenate([p['y'][:,gi] for p in tr])
   sw=np.concatenate([np.full(len(p['y']),1/len(p['y'])) for p in tr])
   m=make_pipeline(StandardScaler(),Ridge(alpha=alpha,solver='lsqr'))
   m.fit(X,y,standardscaler__sample_weight=sw,ridge__sample_weight=sw)
   pred=np.maximum(0,m.predict(vh['features'][feature])); vals.append(float(np.mean(np.abs(pred-vh['y'][:,gi]))))
  scores.append(float(np.mean(vals)))
 best=alphas[int(np.argmin(scores))]
 X=np.vstack([p['features'][feature] for p in train]); y=np.concatenate([p['y'][:,gi] for p in train]); sw=np.concatenate([np.full(len(p['y']),1/len(p['y'])) for p in train])
 m=make_pipeline(StandardScaler(),Ridge(alpha=best,solver='lsqr')); m.fit(X,y,standardscaler__sample_weight=sw,ridge__sample_weight=sw)
 pred=np.maximum(0,m.predict(test['features'][feature])); yy=test['y'][:,gi]
 return {'alpha':best,'inner_macro_MAE':min(scores),'test_MAE':float(np.mean(np.abs(pred-yy))),
  'test_RMSE':float(np.sqrt(np.mean((pred-yy)**2))),'n_test':len(yy)}

rows=[]
for scenario in ['H0','H1']:
 data=build(scenario); dynamic=[v for h,v in data.items() if h!='6B03205']
 for test in dynamic:
  train=[p for p in dynamic if p['heat']!=test['heat']]
  for gi,gas in enumerate(['CO2','CO']):
   for model in ['M0_persistence','M1_current_flow','M2_distributed_lag','M3_history_windows']:
    r=fit_eval(train,test,model,gi); rows.append({'scenario':scenario,'test_heat':test['heat'],'gas':gas,'model':model,**r})
 print('completed LOHO scenario',scenario,flush=True)
df=pd.DataFrame(rows); df.to_csv(tables/'distributed_lag_LOHO_fair.csv',index=False)
macro=df.groupby(['scenario','gas','model']).agg(MAE_macro=('test_MAE','mean'),RMSE_macro=('test_RMSE','mean'),n_heats=('test_heat','nunique')).reset_index()
macro.to_csv(tables/'distributed_lag_LOHO_fair_macro.csv',index=False)
print(macro.to_string(index=False))
