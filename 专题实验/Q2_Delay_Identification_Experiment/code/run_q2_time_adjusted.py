from pathlib import Path
import numpy as np, pandas as pd, math, json, zipfile
from scipy.signal import savgol_filter, correlate
from scipy.optimize import lsq_linear
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

ROOT=Path(__file__).resolve().parents[1]; CURVE_DIR=ROOT/'inputs'; OUT=ROOT/'rerun_output'
HEATS=['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108']; GASES=['CO2','CO']
curves={h:pd.read_csv(CURVE_DIR/f'curves_{h}.csv') for h in HEATS}

def roll_slope(y,dt,win=10):
 n=len(y);k=max(5,int(round(win/dt)));k+=1-(k%2);half=k//2;out=np.full(n,np.nan)
 for i in range(n):
  a=max(0,i-half);b=min(n,i+half+1)
  if b-a<4:continue
  xx=np.arange(a,b)*dt;yy=y[a:b];xc=xx-xx.mean();den=xc@xc
  if den>0:out[i]=max(0,(xc@(yy-yy.mean()))/den)
 return pd.Series(out).interpolate(limit_direction='both').to_numpy()
def deriv(y,dt,sec=5):
 w=max(5,int(round(sec/dt)));w+=1-(w%2);w=min(w,len(y) if len(y)%2 else len(y)-1)
 ys=savgol_filter(y,w,2,mode='interp') if w>=5 else y
 return np.gradient(ys,dt)
def corr_profile(x,y,dt,maxlag=40):
 n=min(len(x),len(y));x=np.asarray(x[:n],float);y=np.asarray(y[:n],float);K=min(int(round(maxlag/dt)),n-5)
 cr=correlate(y,x,mode='full',method='fft')[n-1:n+K]
 px=np.r_[0,np.cumsum(x)];py=np.r_[0,np.cumsum(y)];px2=np.r_[0,np.cumsum(x*x)];py2=np.r_[0,np.cumsum(y*y)]
 k=np.arange(K+1);m=n-k;sx=px[n-k];sy=py[n]-py[k];sx2=px2[n-k];sy2=py2[n]-py[k]
 cov=cr-sx*sy/m;den=np.sqrt(np.maximum(0,sx2-sx*sx/m)*np.maximum(0,sy2-sy*sy/m))
 r=np.divide(cov,den,out=np.full_like(cov,np.nan),where=den>1e-12);return k*dt,np.clip(r,-1,1)

rows=[]
for scenario in ['H0','H1']:
 for h in HEATS:
  d=curves[h];dt=1.0 if scenario=='H1' and h=='6C02057' else .5
  q=roll_slope(d.oxygen_cumulative.to_numpy(float),dt,10)
  for gas in GASES:
   g=d[f'{gas}_flow'].to_numpy(float);x=deriv(q,dt,5);y=deriv(g,dt,5)
   wlen=int(round(120/dt));step=wlen
   for st in range(0,len(d)-wlen+1,step):
    en=st+wlen;xx=x[st:en];yy=y[st:en];qq=q[st:en]
    l,r=corr_profile(xx,yy,dt,40)
    if not np.isfinite(r).any():continue
    i=int(np.nanargmax(np.abs(r)));peak=float(r[i]);
    qmed=float(np.median(qq));qsd=float(np.std(qq));mid=((st+en-1)/2*dt);dur=(len(d)-1)*dt;tf=mid/max(dur,1)
    rows.append({'scenario':scenario,'heat':h,'gas':gas,'window_start_s':st*dt,'window_end_s':(en-1)*dt,
      'window_mid_s':mid,'time_fraction':tf,'flow_median_m3s':qmed,'flow_sd_m3s':qsd,
      'best_lag_s':l[i],'best_corr':peak,'best_abs_corr':abs(peak),'best_sign':int(np.sign(peak)),
      'boundary':bool(i==0 or i==len(r)-1),'informative':bool(qsd>=.15 and abs(peak)>=.25 and not(i==0 or i==len(r)-1))})
win=pd.DataFrame(rows)
# Flow support and bins anchored on H0, excluding fifth outlier.
ref=win[(win.scenario=='H0')&(win.heat!='6C02057')]
lo=float(ref.flow_median_m3s.quantile(.01));hi=float(ref.flow_median_m3s.quantile(.99));cuts=ref.flow_median_m3s.quantile([1/3,2/3]).to_numpy()
win['in_reference_flow_support']=win.flow_median_m3s.between(lo,hi)
win['flow_bin']=pd.cut(win.flow_median_m3s,[-np.inf,cuts[0],cuts[1],np.inf],labels=['low','middle','high'])
win.to_csv(OUT/'tables/nonoverlap_120s_time_adjusted_windows.csv',index=False)

# Slope estimates, both raw within-heat and after removing within-heat linear time trend.
def beta_within(d,xcol,ycol='best_lag_s',group='heat'):
 if len(d)<5 or d[group].nunique()<3:return np.nan
 x=d[xcol].to_numpy(float);y=d[ycol].to_numpy(float);g=d[group].to_numpy()
 xr=x.copy();yr=y.copy()
 for h in np.unique(g):
  m=g==h;xr[m]-=xr[m].mean();yr[m]-=yr[m].mean()
 den=xr@xr
 return float(xr@yr/den) if den>0 else np.nan

def time_residuals(d, poly=1):
 z=d.copy(); z['flow_resid']=np.nan;z['lag_resid']=np.nan
 for h,g in z.groupby('heat'):
  t=g.time_fraction.to_numpy(float)
  X=np.column_stack([np.ones(len(t))]+[t**k for k in range(1,poly+1)])
  for col,out in [('flow_median_m3s','flow_resid'),('best_lag_s','lag_resid')]:
   y=g[col].to_numpy(float)
   if len(y)>X.shape[1]:
    coef=np.linalg.lstsq(X,y,rcond=None)[0];z.loc[g.index,out]=y-X@coef
 return z

slope=[];per=[]
for scenario in ['H0','H1']:
 for gas in GASES:
  z=win[(win.scenario==scenario)&(win.gas==gas)&win.informative&win.in_reference_flow_support].copy()
  if scenario=='H0':z=z[z.heat!='6C02057']
  raw=beta_within(z,'flow_median_m3s')
  tr=time_residuals(z,poly=1)
  adj=beta_within(tr.dropna(subset=['flow_resid','lag_resid']),'flow_resid','lag_resid')
  tr2=time_residuals(z,poly=2)
  adj2=beta_within(tr2.dropna(subset=['flow_resid','lag_resid']),'flow_resid','lag_resid')
  # Heat-cluster bootstrap for time-adjusted beta.
  hs=sorted(z.heat.unique()); groups={h:z[z.heat==h].copy() for h in hs}; rng=np.random.default_rng(20261126+len(slope)); boots=[]
  for _ in range(1000):
   samp=rng.choice(hs,size=len(hs),replace=True); pieces=[]
   for j,h in enumerate(samp):
    g=groups[h].copy();g['boot_heat']=j;pieces.append(g)
   b=pd.concat(pieces,ignore_index=True)
   # residualize time within each bootstrap cluster
   b['flow_resid']=np.nan;b['lag_resid']=np.nan
   for j,g in b.groupby('boot_heat'):
    t=g.time_fraction.to_numpy(float);X=np.column_stack([np.ones(len(t)),t])
    for col,out in [('flow_median_m3s','flow_resid'),('best_lag_s','lag_resid')]:
     y=g[col].to_numpy(float)
     if len(y)>2:
      c=np.linalg.lstsq(X,y,rcond=None)[0];b.loc[g.index,out]=y-X@c
   b=b.dropna(subset=['flow_resid','lag_resid']);xx=b.flow_resid.to_numpy();yy=b.lag_resid.to_numpy();den=xx@xx
   if den>0:boots.append(xx@yy/den)
  ci=np.quantile(boots,[.025,.5,.975]) if boots else [np.nan]*3
  # Per-heat raw and time-adjusted rank correlations.
  for h,g in z.groupby('heat'):
   gt=time_residuals(g,poly=1)
   rho=spearmanr(g.flow_median_m3s,g.best_lag_s).statistic if len(g)>=3 and g.flow_median_m3s.nunique()>1 and g.best_lag_s.nunique()>1 else np.nan
   rho_adj=spearmanr(gt.flow_resid,gt.lag_resid).statistic if len(g)>=3 and gt.flow_resid.nunique()>1 and gt.lag_resid.nunique()>1 else np.nan
   per.append({'scenario':scenario,'gas':gas,'heat':h,'n_windows':len(g),'raw_spearman':rho,'time_adjusted_spearman':rho_adj})
  slope.append({'scenario':scenario,'gas':gas,'n_windows':len(z),'n_heats':len(hs),'raw_within_heat_beta':raw,
   'time_adjusted_beta_linear_time':adj,'time_adjusted_beta_quadratic_time':adj2,
   'time_adjusted_boot_ci_low':ci[0],'time_adjusted_boot_median':ci[1],'time_adjusted_boot_ci_high':ci[2],
   'bootstrap_share_beta_negative':float(np.mean(np.asarray(boots)<0)) if boots else np.nan})

sl=pd.DataFrame(slope);ph=pd.DataFrame(per);sl.to_csv(OUT/'tables/nonoverlap_time_adjusted_flow_delay_slopes.csv',index=False);ph.to_csv(OUT/'tables/nonoverlap_per_heat_time_adjusted_rho.csv',index=False)

# LOHO comparisons: time-only, flow-only, time+flow unconstrained, time+flow with beta_flow <=0.
models=[]
for scenario in ['H0','H1']:
 for gas in GASES:
  z=win[(win.scenario==scenario)&(win.gas==gas)&win.informative&win.in_reference_flow_support].copy()
  if scenario=='H0':z=z[z.heat!='6C02057']
  hs=sorted(z.heat.unique())
  for testh in hs:
   tr=z[z.heat!=testh];te=z[z.heat==testh]
   if len(te)==0 or tr.heat.nunique()<3:continue
   ytr=tr.best_lag_s.to_numpy();yte=te.best_lag_s.to_numpy()
   # baseline
   pred0=np.repeat(np.median(ytr),len(te))
   predt=LinearRegression().fit(tr[['time_fraction']],ytr).predict(te[['time_fraction']])
   predf=LinearRegression().fit(tr[['flow_median_m3s']],ytr).predict(te[['flow_median_m3s']])
   predtf=LinearRegression().fit(tr[['time_fraction','flow_median_m3s']],ytr).predict(te[['time_fraction','flow_median_m3s']])
   Xtr=np.column_stack([np.ones(len(tr)),tr.time_fraction.to_numpy(),tr.flow_median_m3s.to_numpy()])
   Xte=np.column_stack([np.ones(len(te)),te.time_fraction.to_numpy(),te.flow_median_m3s.to_numpy()])
   fit=lsq_linear(Xtr,ytr,bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0.0]))
   predcon=Xte@fit.x
   for name,pred in [('fixed_median',pred0),('time_only',predt),('flow_only',predf),('time_plus_flow_free',predtf),('time_plus_flow_nonincreasing',predcon)]:
    models.append({'scenario':scenario,'gas':gas,'test_heat':testh,'model':name,'n_test_windows':len(te),
     'MAE_s':float(np.mean(np.abs(pred-yte))),'RMSE_s':float(np.sqrt(np.mean((pred-yte)**2)))})
md=pd.DataFrame(models);md.to_csv(OUT/'tables/nonoverlap_loho_time_flow_models.csv',index=False)

# Flow-bin summaries.
binrows=[]
for scenario in ['H0','H1']:
 for gas in GASES:
  z=win[(win.scenario==scenario)&(win.gas==gas)&win.informative]
  for b in ['low','middle','high']:
   q=z[z.flow_bin==b]
   binrows.append({'scenario':scenario,'gas':gas,'flow_bin':b,'n_windows':len(q),'n_heats':q.heat.nunique(),
    'flow_median':q.flow_median_m3s.median() if len(q) else np.nan,'lag_median':q.best_lag_s.median() if len(q) else np.nan,
    'lag_q25':q.best_lag_s.quantile(.25) if len(q) else np.nan,'lag_q75':q.best_lag_s.quantile(.75) if len(q) else np.nan})
pd.DataFrame(binrows).to_csv(OUT/'tables/nonoverlap_flow_bin_lags.csv',index=False)

# Figures: time-adjusted scatter and LOHO comparison.
for scenario in ['H0','H1']:
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
 for ax,gas in zip(axs,GASES):
  z=win[(win.scenario==scenario)&(win.gas==gas)&win.informative]
  for h,g in z.groupby('heat'):
   ax.scatter(g.flow_median_m3s,g.best_lag_s,s=35,alpha=.7,label=h)
  ax.set_xlabel('Median oxygen flow (m³/s)');ax.set_ylabel('Local lag (s)');ax.set_title(gas);ax.grid(True,alpha=.25);ax.legend(fontsize=7,ncol=2,frameon=False)
 fig.suptitle(f'Non-overlapping 120-s windows: flow vs lag — {scenario}',fontsize=14)
 fig.savefig(OUT/f'figures/{scenario}_nonoverlap_flow_lag.png',dpi=200,bbox_inches='tight');plt.close(fig)

# Update report addendum.
add=['','## 9. 新增：非重叠窗口与吹炼进程控制','',
'前一版 120 s 窗口采用 60 s 步长，存在 50% 重叠，且流量和延迟都可能随吹炼进程变化。本轮新增 120 s 非重叠窗口，并显式控制每炉内的时间进程。',
'- 原始炉内流量—延迟斜率：见 `nonoverlap_time_adjusted_flow_delay_slopes.csv` 的 `raw_within_heat_beta`。',
'- 时间控制斜率：先在每炉内分别用时间比例线性/二次趋势解释流量和延迟，再对残差回归；线性时间控制为主，二次为敏感性。',
'- 留一炉预测比较固定中位数、仅时间、仅流量、时间+流量自由模型、时间+流量且流量系数非正的约束模型。若加入流量未优于仅时间，不能说流量提供了独立解释。', '',
'### 新增结果',
'']
for r in slope:
 add.append(f"- {r['scenario']} / {r['gas']}：原始炉内斜率 {r['raw_within_heat_beta']:.3f}；控制线性时间趋势后斜率 {r['time_adjusted_beta_linear_time']:.3f}，炉级 bootstrap 95%区间 [{r['time_adjusted_boot_ci_low']:.3f}, {r['time_adjusted_boot_ci_high']:.3f}]；控制二次时间趋势后斜率 {r['time_adjusted_beta_quadratic_time']:.3f}。")
add.append('')
if len(md):
 add.append('### 新增 LOHO 误差（MAE，秒）')
 for idx,v in md.groupby(['scenario','gas','model']).MAE_s.mean().items():add.append(f'- {idx[0]} / {idx[1]} / {idx[2]}：{v:.2f}')
add += ['', '### 解释原则',
'若原始斜率为负、但控制时间趋势后斜率明显衰减或区间跨零，说明先前的负斜率很可能混有吹炼进程效应；即使控制后仍为负，也要看时间+流量模型是否在留一炉测试中优于仅时间模型。',
'本附加实验仍是窗口级探索性识别，不构成物理因果证明。']
with open(OUT/'README_Q2_Delay_Experiment.md','a',encoding='utf-8') as f:f.write('\n'.join(add)+'\n')

# Rebuild ZIP.
zip_path=ROOT.parent/'Q2_Delay_Identification_Experiment_rerun.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():zf.write(p,arcname=str(p.relative_to(OUT)))
print('UPDATED ZIP',zip_path,zip_path.stat().st_size)
print('\nTIME ADJUSTED SLOPES')
print(sl.round(3).to_string(index=False))
print('\nPER HEAT RHO')
print(ph.round(3).to_string(index=False))
print('\nLOHO MODEL MAE')
print(md.groupby(['scenario','gas','model']).MAE_s.mean().round(2).to_string())
print('\nNONOVERLAP BINS')
print(pd.DataFrame(binrows).round(2).to_string(index=False))
