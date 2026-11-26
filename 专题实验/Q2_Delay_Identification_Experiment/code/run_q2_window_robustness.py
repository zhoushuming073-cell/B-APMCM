from pathlib import Path
import numpy as np,pandas as pd,math,zipfile
from scipy.signal import savgol_filter,correlate
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1];CURVE_DIR=ROOT/'inputs';OUT=ROOT/'rerun_output'
HEATS=['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108'];GASES=['CO2','CO']
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
def corrprof(x,y,dt,maxlag):
 n=min(len(x),len(y));x=np.asarray(x[:n],float);y=np.asarray(y[:n],float);K=min(int(round(maxlag/dt)),n-5)
 cr=correlate(y,x,mode='full',method='fft')[n-1:n+K]
 px=np.r_[0,np.cumsum(x)];py=np.r_[0,np.cumsum(y)];px2=np.r_[0,np.cumsum(x*x)];py2=np.r_[0,np.cumsum(y*y)]
 k=np.arange(K+1);m=n-k;sx=px[n-k];sy=py[n]-py[k];sx2=px2[n-k];sy2=py2[n]-py2[k]
 cov=cr-sx*sy/m;den=np.sqrt(np.maximum(0,sx2-sx*sx/m)*np.maximum(0,sy2-sy*sy/m))
 r=np.divide(cov,den,out=np.full_like(cov,np.nan),where=den>1e-12)
 return k*dt,np.clip(r,-1,1)
def shift_p(x,y,dt,maxlag,seed):
 l,r=corrprof(x,y,dt,maxlag);obs=np.nanmax(np.abs(r));n=len(x);k=int(round(maxlag/dt));margin=k+int(round(5/dt));lo=margin;hi=n-margin
 if hi<=lo:return np.nan,obs
 shifts=np.arange(lo,hi)
 if len(shifts)>80:
  rng=np.random.default_rng(seed);shifts=rng.choice(shifts,80,replace=False)
 null=[]
 for s in shifts:
  _,rr=corrprof(x,np.roll(y,int(s)),dt,maxlag)
  null.append(np.nanmax(np.abs(rr)))
 return (1+np.sum(np.asarray(null)>=obs))/(1+len(null)),obs

def resid(g,col,poly=1):
 t=g.time_fraction.to_numpy(float);X=np.column_stack([np.ones(len(t))]+[t**k for k in range(1,poly+1)]);y=g[col].to_numpy(float)
 if len(y)<=X.shape[1]:return np.full(len(y),np.nan)
 return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def beta(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float);d=x@x;return x@y/d if d>0 else np.nan

allrows=[]
for sc in ['H0','H1']:
 for h in HEATS:
  d=curves[h];dt=1.0 if sc=='H1' and h=='6C02057' else .5;q=roll_slope(d.oxygen_cumulative.to_numpy(float),dt);dur=(len(d)-1)*dt
  for gas in GASES:
   x=deriv(q,dt);y=deriv(d[f'{gas}_flow'].to_numpy(float),dt)
   for win_s,maxlag in [(90.,30.),(180.,60.)]:
    nwin=int(round(win_s/dt))
    for st in range(0,len(d)-nwin+1,nwin):
     en=st+nwin;xx=x[st:en];yy=y[st:en];qq=q[st:en]
     ll,rr=corrprof(xx,yy,dt,maxlag)
     if not np.isfinite(rr).any():continue
     i=int(np.nanargmax(np.abs(rr)));pk=float(rr[i]);mid=((st+en-1)/2)*dt
     p,obs=shift_p(xx,yy,dt,maxlag,20261126+st+len(allrows))
     allrows.append({'scenario':sc,'window_s':win_s,'maxlag_s':maxlag,'heat':h,'gas':gas,
      'window_start_s':st*dt,'window_end_s':(en-1)*dt,'window_mid_s':mid,'time_fraction':mid/max(dur,1),
      'flow_median_m3s':float(np.median(qq)),'flow_sd_m3s':float(np.std(qq)),
      'lag_s':float(ll[i]),'peak_r':pk,'peak_abs_r':abs(pk),'sign':int(np.sign(pk)),
      'boundary':bool(i==0 or i==len(rr)-1),'circular_p':p})
win=pd.DataFrame(allrows)
# BH by scenario, window size and gas
win['BH_q']=np.nan
for _,g in win.groupby(['scenario','window_s','gas']):
 ix=g.index;valid=g.circular_p.notna()
 if valid.any():win.loc[ix[valid],'BH_q']=multipletests(g.loc[valid,'circular_p'],method='fdr_bh')[1]
win['threshold_informative']=(win.flow_sd_m3s>=.15)&(win.peak_abs_r>=.25)&(~win.boundary)
win['statistically_supported']=(win.threshold_informative)&(win.BH_q<=.10)
# Common flow support and bins anchored to H0 90s windows, excluding fifth.
ref=win[(win.scenario=='H0')&(win.window_s==90)&(win.heat!='6C02057')]
lo,hi=ref.flow_median_m3s.quantile([.01,.99]);cuts=ref.flow_median_m3s.quantile([1/3,2/3]).to_numpy()
win['in_reference_support']=win.flow_median_m3s.between(lo,hi)
win['flow_bin']=pd.cut(win.flow_median_m3s,[-np.inf,cuts[0],cuts[1],np.inf],labels=['low','middle','high'])
win.to_csv(OUT/'tables/nonoverlap_window_size_robustness.csv',index=False)

# Flow-delay slopes for 90s windows: thresholded vs q-supported; time adjusted.
sl=[];model=[];per=[]
for sc in ['H0','H1']:
 for gas in GASES:
  for selname,sel in [('threshold',win.threshold_informative),('BH_supported',win.statistically_supported)]:
   z=win[(win.scenario==sc)&(win.window_s==90)&(win.gas==gas)&sel&win.in_reference_support].copy()
   if sc=='H0':z=z[z.heat!='6C02057']
   if len(z)<5 or z.heat.nunique()<3:
    sl.append({'scenario':sc,'gas':gas,'selection':selname,'n_windows':len(z),'n_heats':z.heat.nunique(),'raw_beta':np.nan,'time_adjusted_beta':np.nan,'boot_low':np.nan,'boot_high':np.nan})
    continue
   # raw within-heat beta
   xx=z.flow_median_m3s.to_numpy(float);yy=z.lag_s.to_numpy(float);gg=z.heat.to_numpy();xr=xx.copy();yr=yy.copy()
   for h in np.unique(gg):m=gg==h;xr[m]-=xr[m].mean();yr[m]-=yr[m].mean()
   raw=beta(xr,yr)
   z['xr']=np.nan;z['yr']=np.nan
   arrays={}
   for h,g in z.groupby('heat'):
    a=resid(g,'flow_median_m3s',1);b=resid(g,'lag_s',1);z.loc[g.index,'xr']=a;z.loc[g.index,'yr']=b;arrays[h]=(a,b)
   adj=beta(z.xr.fillna(0),z.yr.fillna(0))
   hs=list(arrays);rng=np.random.default_rng(20261126+len(sl));boots=[]
   for _ in range(1000):
    samp=rng.choice(hs,size=len(hs),replace=True);a=[];b=[]
    for h in samp:
     xh,yh=arrays[h];ok=np.isfinite(xh)&np.isfinite(yh);a.extend(xh[ok]);b.extend(yh[ok])
    boots.append(beta(a,b))
   ci=np.nanquantile(boots,[.025,.975])
   sl.append({'scenario':sc,'gas':gas,'selection':selname,'n_windows':len(z),'n_heats':len(hs),'raw_beta':raw,'time_adjusted_beta':adj,'boot_low':ci[0],'boot_high':ci[1],'boot_share_negative':np.nanmean(np.asarray(boots)<0)})
   for h,g in z.groupby('heat'):
    rho=spearmanr(g.flow_median_m3s,g.lag_s).statistic if len(g)>=3 and g.flow_median_m3s.nunique()>1 and g.lag_s.nunique()>1 else np.nan
    per.append({'scenario':sc,'gas':gas,'selection':selname,'heat':h,'n_windows':len(g),'rho':rho})
  # LOHO with threshold-informative windows (not BH-filtered), include time baseline.
  z=win[(win.scenario==sc)&(win.window_s==90)&(win.gas==gas)&win.threshold_informative&win.in_reference_support].copy()
  if sc=='H0':z=z[z.heat!='6C02057']
  hs=sorted(z.heat.unique())
  for testh in hs:
   tr=z[z.heat!=testh];te=z[z.heat==testh]
   if len(te)==0 or tr.heat.nunique()<3:continue
   yt=tr.lag_s.to_numpy();yv=te.lag_s.to_numpy()
   preds={'fixed_median':np.repeat(np.median(yt),len(te)),
    'time_only':LinearRegression().fit(tr[['time_fraction']],yt).predict(te[['time_fraction']]),
    'flow_only':LinearRegression().fit(tr[['flow_median_m3s']],yt).predict(te[['flow_median_m3s']]),
    'time_plus_flow':LinearRegression().fit(tr[['time_fraction','flow_median_m3s']],yt).predict(te[['time_fraction','flow_median_m3s']])}
   Xtr=np.column_stack([np.ones(len(tr)),tr.time_fraction,tr.flow_median_m3s]);Xte=np.column_stack([np.ones(len(te)),te.time_fraction,te.flow_median_m3s])
   preds['time_plus_flow_noninc']=Xte@lsq_linear(Xtr,yt,bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,0])).x
   for name,p in preds.items():model.append({'scenario':sc,'gas':gas,'test_heat':testh,'model':name,'n_test':len(te),'MAE_s':np.mean(np.abs(p-yv))})

sl=pd.DataFrame(sl);per=pd.DataFrame(per);md=pd.DataFrame(model)
sl.to_csv(OUT/'tables/90s_window_flow_delay_slopes.csv',index=False);per.to_csv(OUT/'tables/90s_window_per_heat_rho.csv',index=False);md.to_csv(OUT/'tables/90s_window_loho_models.csv',index=False)
# Summaries for lag p values and counts
cnt=win[win.window_s==90].groupby(['scenario','gas']).agg(windows=('heat','size'),threshold_informative=('threshold_informative','sum'),BH_supported=('statistically_supported','sum'),heats=('heat','nunique'),median_peak_abs_r=('peak_abs_r','median'),median_lag=('lag_s','median')).reset_index()
cnt.to_csv(OUT/'tables/90s_window_support_counts.csv',index=False)
# append interpretation report
with open(OUT/'README_Q2_Delay_Experiment.md','a',encoding='utf-8') as f:
 f.write('\n## 11. 新增：90 秒非重叠窗口与窗口级负对照\n\n')
 f.write('将窗口缩短至 90 秒、非重叠，滞后搜索范围缩短至 0–30 秒；每个窗口用循环移位保持自相关的近似负对照，并在“情景×窗口长度×气体”内进行 BH 校正。此检验的有效移位数量有限，仍是局部近似负对照。\n')
 f.write('\n窗口支持计数：\n')
 for r in cnt.to_dict('records'):f.write(f"- {r['scenario']} / {r['gas']}：{r['windows']} 个窗口；阈值信息性 {r['threshold_informative']} 个；BH q≤0.10 {r['BH_supported']} 个；中位 |r| {r['median_peak_abs_r']:.3f}；中位 lag {r['median_lag']:.1f}s。\n")
 f.write('\n流量—延迟斜率结果：\n')
 for r in sl.to_dict('records'):
  if pd.notna(r.get('time_adjusted_beta')):f.write(f"- {r['scenario']} / {r['gas']} / {r['selection']}：原始斜率 {r['raw_beta']:.3f}；控制时间趋势斜率 {r['time_adjusted_beta']:.3f}，炉级 bootstrap 95%区间 [{r['boot_low']:.3f}, {r['boot_high']:.3f}]。\n")
  else:f.write(f"- {r['scenario']} / {r['gas']} / {r['selection']}：仅 {r['n_windows']} 个窗口、{r['n_heats']} 炉，样本不足，不估斜率。\n")
 if len(md):
  f.write('\n90 秒窗口留一炉 MAE：\n')
  for idx,v in md.groupby(['scenario','gas','model']).MAE_s.mean().items():f.write(f'- {idx[0]} / {idx[1]} / {idx[2]}：{v:.2f}s\n')
 f.write('\n注意：每炉仅约 4–9 个非重叠窗口，时间趋势控制仍受小样本限制。只有 BH 支持的窗口足够多且跨炉复现时，才可将流量—延迟关系升级为主结论。\n')
# rebuild zip
zp=ROOT.parent/'Q2_Delay_Identification_Experiment_rerun.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():zf.write(p,arcname=str(p.relative_to(OUT)))
print('ZIP',zp,zp.stat().st_size)
print('\nCOUNTS\n',cnt.round(3).to_string(index=False))
print('\nSLOPES\n',sl.round(3).to_string(index=False))
print('\nLOHO\n',md.groupby(['scenario','gas','model']).MAE_s.mean().round(2).to_string())
print('\nWINDOW EXAMPLES\n',win[win.window_s==90].head(20).round(3).to_string(index=False))
