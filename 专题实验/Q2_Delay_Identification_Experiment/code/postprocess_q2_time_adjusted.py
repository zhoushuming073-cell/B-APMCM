from pathlib import Path
import numpy as np,pandas as pd,json,zipfile
from sklearn.linear_model import LinearRegression
from scipy.optimize import lsq_linear
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'rerun_output'
win=pd.read_csv(OUT/'tables/nonoverlap_120s_time_adjusted_windows.csv');GASES=['CO2','CO']

def residualize(g,col,poly=1):
 t=g.time_fraction.to_numpy(float);X=np.column_stack([np.ones(len(t))]+[t**k for k in range(1,poly+1)])
 y=g[col].to_numpy(float)
 if len(y)<=X.shape[1]:return np.full(len(y),np.nan)
 c=np.linalg.lstsq(X,y,rcond=None)[0];return y-X@c

def slope_xy(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float);den=x@x
 return float(x@y/den) if den>0 else np.nan

sl=[];per=[];models=[];binrows=[]
for sc in ['H0','H1']:
 for gas in GASES:
  z=win[(win.scenario==sc)&(win.gas==gas)&win.informative&win.in_reference_flow_support].copy()
  if sc=='H0':z=z[z.heat!='6C02057']
  # within-heat demeaned raw slope
  x=z.flow_median_m3s.to_numpy(float);y=z.best_lag_s.to_numpy(float);grp=z.heat.to_numpy();xr=x.copy();yr=y.copy()
  for h in np.unique(grp):
   m=grp==h;xr[m]-=xr[m].mean();yr[m]-=yr[m].mean()
  raw=slope_xy(xr,yr)
  # Time residualization within each heat, linear and quadratic.
  z['flow_r1']=np.nan;z['lag_r1']=np.nan;z['flow_r2']=np.nan;z['lag_r2']=np.nan
  arrays={}
  for h,g in z.groupby('heat'):
   f1=residualize(g,'flow_median_m3s',1);l1=residualize(g,'best_lag_s',1)
   f2=residualize(g,'flow_median_m3s',2);l2=residualize(g,'best_lag_s',2)
   z.loc[g.index,'flow_r1']=f1;z.loc[g.index,'lag_r1']=l1;z.loc[g.index,'flow_r2']=f2;z.loc[g.index,'lag_r2']=l2
   arrays[h]=(f1,l1,f2,l2)
  adj1=slope_xy(z.flow_r1.fillna(0),z.lag_r1.fillna(0));adj2=slope_xy(z.flow_r2.fillna(0),z.lag_r2.fillna(0))
  # Heat-cluster bootstrap from precomputed within-heat residuals.
  hs=list(arrays);rng=np.random.default_rng(20261126+len(sl));boots=[]
  for _ in range(1000):
   sample=rng.choice(hs,size=len(hs),replace=True);xs=[];ys=[]
   for h in sample:
    a,b,_,_=arrays[h];ok=np.isfinite(a)&np.isfinite(b);xs.extend(a[ok]);ys.extend(b[ok])
   boots.append(slope_xy(xs,ys))
  boots=np.asarray(boots);ci=np.nanquantile(boots,[.025,.5,.975])
  sl.append({'scenario':sc,'gas':gas,'n_windows':len(z),'n_heats':len(hs),'raw_within_heat_beta':raw,
    'time_adjusted_beta_linear_time':adj1,'time_adjusted_beta_quadratic_time':adj2,
    'time_adjusted_boot_ci_low':ci[0],'time_adjusted_boot_median':ci[1],'time_adjusted_boot_ci_high':ci[2],
    'bootstrap_share_beta_negative':np.nanmean(boots<0)})
  for h,g in z.groupby('heat'):
   rho=spearmanr(g.flow_median_m3s,g.best_lag_s).statistic if len(g)>=3 and g.flow_median_m3s.nunique()>1 and g.best_lag_s.nunique()>1 else np.nan
   gr=g.copy();f=residualize(gr,'flow_median_m3s',1);l=residualize(gr,'best_lag_s',1)
   ra=spearmanr(f,l).statistic if len(g)>=3 and np.nanstd(f)>1e-10 and np.nanstd(l)>1e-10 else np.nan
   per.append({'scenario':sc,'gas':gas,'heat':h,'n_windows':len(g),'raw_spearman':rho,'time_adjusted_spearman':ra})
  # LOHO: fixed, time only, flow only, time+flow, constrained flow coefficient <=0.
  for testh in hs:
   tr=z[z.heat!=testh];te=z[z.heat==testh]
   if len(te)==0 or tr.heat.nunique()<3:continue
   yt=tr.best_lag_s.to_numpy(float);yv=te.best_lag_s.to_numpy(float)
   p0=np.repeat(np.median(yt),len(te))
   pt=LinearRegression().fit(tr[['time_fraction']],yt).predict(te[['time_fraction']])
   pf=LinearRegression().fit(tr[['flow_median_m3s']],yt).predict(te[['flow_median_m3s']])
   ptf=LinearRegression().fit(tr[['time_fraction','flow_median_m3s']],yt).predict(te[['time_fraction','flow_median_m3s']])
   Xtr=np.column_stack([np.ones(len(tr)),tr.time_fraction.to_numpy(),tr.flow_median_m3s.to_numpy()])
   Xte=np.column_stack([np.ones(len(te)),te.time_fraction.to_numpy(),te.flow_median_m3s.to_numpy()])
   fit=lsq_linear(Xtr,yt,bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0]))
   pc=Xte@fit.x
   for name,pred in [('fixed_median',p0),('time_only',pt),('flow_only',pf),('time_plus_flow_free',ptf),('time_plus_flow_nonincreasing',pc)]:
    models.append({'scenario':sc,'gas':gas,'test_heat':testh,'model':name,'n_test_windows':len(te),'MAE_s':np.mean(np.abs(pred-yv)),'RMSE_s':np.sqrt(np.mean((pred-yv)**2))})
  # flow bins
  for b in ['low','middle','high']:
   a=win[(win.scenario==sc)&(win.gas==gas)&win.informative&(win.flow_bin==b)]
   binrows.append({'scenario':sc,'gas':gas,'flow_bin':b,'n_windows':len(a),'n_heats':a.heat.nunique(),'flow_median':a.flow_median_m3s.median() if len(a) else np.nan,'lag_median':a.best_lag_s.median() if len(a) else np.nan,'lag_q25':a.best_lag_s.quantile(.25) if len(a) else np.nan,'lag_q75':a.best_lag_s.quantile(.75) if len(a) else np.nan})

sl=pd.DataFrame(sl);per=pd.DataFrame(per);md=pd.DataFrame(models);bins=pd.DataFrame(binrows)
sl.to_csv(OUT/'tables/nonoverlap_time_adjusted_flow_delay_slopes.csv',index=False);per.to_csv(OUT/'tables/nonoverlap_per_heat_time_adjusted_rho.csv',index=False);md.to_csv(OUT/'tables/nonoverlap_loho_time_flow_models.csv',index=False);bins.to_csv(OUT/'tables/nonoverlap_flow_bin_lags.csv',index=False)

# Add final interpretation to report.
with open(OUT/'README_Q2_Delay_Experiment.md','a',encoding='utf-8') as f:
 f.write('\n## 10. 非重叠窗口时间控制实验结果（最终补充）\n\n')
 for r in sl.to_dict('records'):
  f.write(f"- {r['scenario']} / {r['gas']}：原始炉内斜率 {r['raw_within_heat_beta']:.3f}；控制线性时间趋势后斜率 {r['time_adjusted_beta_linear_time']:.3f}，炉级 bootstrap 95%区间 [{r['time_adjusted_boot_ci_low']:.3f}, {r['time_adjusted_boot_ci_high']:.3f}]；控制二次时间趋势后斜率 {r['time_adjusted_beta_quadratic_time']:.3f}。\n")
 f.write('\n留一炉预测 MAE（秒）：\n')
 for idx,v in md.groupby(['scenario','gas','model']).MAE_s.mean().items():f.write(f'- {idx[0]} / {idx[1]} / {idx[2]}：{v:.2f}\n')
 f.write('\n解释：如果加入流量变量未能在留一炉验证中优于仅时间进程模型，则不能认为流量对延迟具有独立的跨炉预测价值。时间调整斜率仍是小样本窗口级关联，不构成因果效应。\n')

# Rebuild ZIP
zp=ROOT.parent/'Q2_Delay_Identification_Experiment_rerun.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():zf.write(p,arcname=str(p.relative_to(OUT)))
print('ZIP',zp,zp.stat().st_size)
print('\nSLOPES\n',sl.round(3).to_string(index=False))
print('\nPER HEAT\n',per.round(3).to_string(index=False))
print('\nLOHO\n',md.groupby(['scenario','gas','model']).MAE_s.mean().round(2).to_string())
print('\nBINS\n',bins.round(2).to_string(index=False))
