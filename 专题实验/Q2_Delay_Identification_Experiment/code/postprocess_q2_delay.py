from pathlib import Path
import numpy as np, pandas as pd, json, math, shutil, zipfile, hashlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, correlate
from sklearn.linear_model import LinearRegression
from sklearn.isotonic import IsotonicRegression
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'rerun_output'
HEATS=['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108']
GASES=['CO2','CO']

# Restore scenario label lost by the groupby/apply serialization in the interrupted first pass.
full=pd.read_csv(OUT/'tables/full_series_lag_diagnostics.csv')
if 'scenario' not in full:
    full['scenario']=['H0']*20+['H1']*20
for sc,g in full.groupby('scenario'):
    ix=g.index
    valid=full.loc[ix,'circular_shift_p'].notna()
    full.loc[ix,'BH_q']=np.nan
    if valid.sum(): full.loc[ix[valid], 'BH_q']=multipletests(full.loc[ix[valid],'circular_shift_p'].values,method='fdr_bh')[1]
full.to_csv(OUT/'tables/full_series_lag_diagnostics.csv',index=False)

windows=pd.read_csv(OUT/'tables/local_window_lag_diagnostics.csv')

# Recompute summaries and slopes with heat-cluster bootstrap; windows are not treated as independent heats.
def within_slope(d):
    d=d.dropna(subset=['flow_median_m3s','best_lag_s'])
    if len(d)<5 or d.heat.nunique()<3:return np.nan
    x=d.flow_median_m3s.to_numpy(float); y=d.best_lag_s.to_numpy(float)
    hs=d.heat.to_numpy()
    xc=x.copy(); yc=y.copy()
    for h in np.unique(hs):
        m=hs==h; xc[m]-=xc[m].mean(); yc[m]-=yc[m].mean()
    den=np.dot(xc,xc)
    return float(np.dot(xc,yc)/den) if den>0 else np.nan

slope_rows=[]; per_rows=[]; model_rows=[]
for scenario in ['H0','H1']:
    for gas in GASES:
        z=windows[(windows.scenario==scenario)&(windows.gas==gas)&windows.informative&windows.in_reference_flow_support].copy()
        if scenario=='H0': z=z[z.heat!='6C02057']
        hs=sorted(z.heat.unique())
        for h,g in z.groupby('heat'):
            rho=spearmanr(g.flow_median_m3s,g.best_lag_s).statistic if len(g)>=3 and g.flow_median_m3s.nunique()>1 and g.best_lag_s.nunique()>1 else np.nan
            slope=np.polyfit(g.flow_median_m3s,g.best_lag_s,1)[0] if len(g)>=3 and g.flow_median_m3s.nunique()>1 else np.nan
            per_rows.append({'scenario':scenario,'gas':gas,'heat':h,'n_windows':len(g),'spearman_rho':rho,'linear_slope_s_per_m3s':slope})
        beta=within_slope(z)
        # Cluster bootstrap with pre-grouped arrays; 1000 resamples of heats.
        rng=np.random.default_rng(20261126+len(slope_rows))
        group_arrays={h:(z[z.heat==h].flow_median_m3s.to_numpy(float),z[z.heat==h].best_lag_s.to_numpy(float)) for h in hs}
        boots=[]
        if len(hs)>=3:
            for _ in range(1000):
                sample=rng.choice(hs,size=len(hs),replace=True)
                xs=[];ys=[];gs=[]
                for j,h in enumerate(sample):
                    x,y=group_arrays[h];xs.extend(x);ys.extend(y);gs.extend([j]*len(x))
                x=np.asarray(xs);y=np.asarray(ys);g=np.asarray(gs)
                xc=x.copy();yc=y.copy()
                for j in np.unique(g):
                    m=g==j;xc[m]-=xc[m].mean();yc[m]-=yc[m].mean()
                den=np.dot(xc,xc)
                if den>0:boots.append(np.dot(xc,yc)/den)
        ci=np.quantile(boots,[.025,.5,.975]) if boots else [np.nan]*3
        loo=[]
        for h in hs: loo.append(within_slope(z[z.heat!=h]))
        per=pd.DataFrame([r for r in per_rows if r['scenario']==scenario and r['gas']==gas])
        slope_rows.append({'scenario':scenario,'gas':gas,'n_informative_windows':len(z),'n_heats':len(hs),
            'within_heat_slope_s_per_m3s':beta,'cluster_boot_ci_low':ci[0],'cluster_boot_median':ci[1],
            'cluster_boot_ci_high':ci[2],'bootstrap_share_beta_negative':float(np.mean(np.asarray(boots)<0)) if boots else np.nan,
            'loo_negative_share':float(np.mean(np.asarray(loo)<0)) if loo else np.nan,
            'median_per_heat_spearman':per.spearman_rho.median() if len(per) else np.nan,
            'heats_with_negative_rho':int((per.spearman_rho<0).sum()) if len(per) else 0,
            'heats_with_positive_rho':int((per.spearman_rho>0).sum()) if len(per) else 0})
        for testh in hs:
            tr=z[z.heat!=testh];te=z[z.heat==testh]
            if len(te)==0 or tr.heat.nunique()<3:continue
            xtr=tr.flow_median_m3s.to_numpy().reshape(-1,1);ytr=tr.best_lag_s.to_numpy()
            xte=te.flow_median_m3s.to_numpy().reshape(-1,1);yte=te.best_lag_s.to_numpy()
            pred_fixed=np.repeat(np.median(ytr),len(te))
            pred_lin=LinearRegression().fit(xtr,ytr).predict(xte)
            pred_iso=IsotonicRegression(increasing=False,out_of_bounds='clip').fit(xtr.ravel(),ytr).predict(xte.ravel())
            for m,pred in [('fixed_median',pred_fixed),('free_linear',pred_lin),('nonincreasing_isotonic',pred_iso)]:
                model_rows.append({'scenario':scenario,'gas':gas,'test_heat':testh,'model':m,'n_test_windows':len(te),
                    'test_MAE_s':np.mean(np.abs(pred-yte)),'test_RMSE_s':np.sqrt(np.mean((pred-yte)**2))})

slopes=pd.DataFrame(slope_rows); per=pd.DataFrame(per_rows); models=pd.DataFrame(model_rows)
slopes.to_csv(OUT/'tables/within_heat_flow_delay_slope.csv',index=False)
per.to_csv(OUT/'tables/per_heat_flow_delay_slopes.csv',index=False)
models.to_csv(OUT/'tables/flow_delay_model_loho.csv',index=False)

# Flow-bin summaries.
bin_rows=[]
for scenario in ['H0','H1']:
 for gas in GASES:
  s=windows[(windows.scenario==scenario)&(windows.gas==gas)]
  for subset,ss in [('all',s),('informative',s[s.informative])]:
   for b in ['low','middle','high']:
    q=ss[ss.flow_bin==b]
    bin_rows.append({'scenario':scenario,'gas':gas,'window_set':subset,'flow_bin':b,'n_windows':len(q),'n_heats':q.heat.nunique(),
      'flow_median_m3s':q.flow_median_m3s.median() if len(q) else np.nan,
      'lag_median_s':q.best_lag_s.median() if len(q) else np.nan,
      'lag_q25_s':q.best_lag_s.quantile(.25) if len(q) else np.nan,
      'lag_q75_s':q.best_lag_s.quantile(.75) if len(q) else np.nan,
      'median_abs_corr':q.best_corr.abs().median() if len(q) else np.nan,
      'negative_corr_share':(q.best_sign<0).mean() if len(q) else np.nan})
bins=pd.DataFrame(bin_rows); bins.to_csv(OUT/'tables/flow_bin_lag_summary.csv',index=False)

# Spectral brightness proxy lag analysis at 1 Hz, with 5/9/15 s smoothing.
specfile=ROOT/'inputs/harmonized_spectra.npz'
spec_rows=[]
if specfile.exists():
 zsp=np.load(specfile); grid=zsp['grid'].astype(float); edges=np.linspace(345,1015,33)
 bidx=[np.where((grid>=a)&(grid<b))[0] for a,b in zip(edges[:-1],edges[1:])]
 def fastcorr(x,y,dt,maxlag=60):
  n=min(len(x),len(y));x=np.asarray(x[:n],float);y=np.asarray(y[:n],float)
  K=min(int(round(maxlag/dt)),n-5); cr=correlate(y,x,mode='full',method='fft')[n-1:n+K]
  px=np.r_[0,np.cumsum(x)];py=np.r_[0,np.cumsum(y)];px2=np.r_[0,np.cumsum(x*x)];py2=np.r_[0,np.cumsum(y*y)]
  k=np.arange(K+1);m=n-k;sx=px[n-k];sy=py[n]-py[k];sx2=px2[n-k];sy2=py2[n]-py2[k]
  cov=cr-sx*sy/m;den=np.sqrt(np.maximum(0,sx2-sx*sx/m)*np.maximum(0,sy2-sy*sy/m))
  r=np.divide(cov,den,out=np.full_like(cov,np.nan),where=den>1e-12);return k*dt,np.clip(r,-1,1)
 def sm_deriv(y,dt,sec):
  w=max(5,int(round(sec/dt)));w+=1-(w%2)
  w=min(w,len(y) if len(y)%2 else len(y)-1)
  ys=savgol_filter(y,w,2,mode='interp') if w>=5 else y
  return np.gradient(ys,dt)
 for scenario in ['H0','H1']:
  for h in HEATS:
   X=zsp[h].astype(np.float32); bm=np.column_stack([X[:,ix].mean(axis=1) for ix in bidx])
   bright=np.log(np.maximum(bm.mean(axis=1),1e-6))
   dt0=1.0 if scenario=='H1' and h=='6C02057' else .5
   stride=max(1,int(round(1.0/dt0)))
   bright=bright[::stride];dt=dt0*stride
   for gas in GASES:
    gv=pd.read_csv(ROOT/'inputs'/f'curves_{h}.csv')[f'{gas}_flow'].to_numpy(float)[::stride]
    for sw in [5.,9.,15.]:
     x=sm_deriv(bright,dt,sw);y=sm_deriv(gv,dt,sw)
     trim=min(int(round(20/dt)),max(0,(len(x)-100)//4));x=x[trim:-trim] if trim else x;y=y[trim:-trim] if trim else y
     lag,r=fastcorr(x,y,dt);ok=np.isfinite(r)
     if ok.any():
      i=int(np.nanargmax(np.abs(r)));j=int(np.nanargmax(r));k=int(np.nanargmin(r))
      spec_rows.append({'scenario':scenario,'heat':h,'gas':gas,'smooth_s':sw,'best_abs_lag_s':lag[i],
       'best_abs_r':r[i],'best_abs_sign':int(np.sign(r[i])),'best_positive_lag_s':lag[j],
       'best_positive_r':r[j],'best_negative_lag_s':lag[k],'best_negative_r':r[k],
       'boundary':bool(i==0 or i==len(r)-1)})
spec=pd.DataFrame(spec_rows)
if len(spec): spec.to_csv(OUT/'tables/spectral_brightness_lag_diagnostics.csv',index=False)

# Quality inventory and source fingerprints.
qrows=[]
for h in HEATS:
 d=pd.read_csv(ROOT/'inputs'/f'curves_{h}.csv')
 for gas in GASES:
  a=d[f'{gas}_flow'].to_numpy(float)
  qrows.append({'heat':h,'gas':gas,'n_rows':len(a),'zero_share':np.mean(a==0),'negative_count':int((a<0).sum()),
   'min':a.min(),'median':np.median(a),'max':a.max(),'sd':a.std(),'unique_values':pd.Series(a).nunique()})
pd.DataFrame(qrows).to_csv(OUT/'tables/gas_signal_quality.csv',index=False)
fp=[]
for p in sorted((OUT/'inputs').glob('*')):
 if p.is_file(): fp.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
# Include source spectrum fingerprint, without bundling its 80MB source file.
if specfile.exists(): fp.append({'file':'source_harmonized_spectra.npz','bytes':specfile.stat().st_size,'sha256':hashlib.sha256(specfile.read_bytes()).hexdigest()})
pd.DataFrame(fp).to_csv(OUT/'tables/input_sha256.csv',index=False)

# Figures.
# Per-heat full lag profiles (recompute from source curves quickly).
def roll_slope(y,dt,win=10):
 n=len(y);k=max(5,int(round(win/dt)));k+=1-(k%2);half=k//2;out=np.full(n,np.nan)
 for i in range(n):
  a=max(0,i-half);b=min(n,i+half+1)
  if b-a<4:continue
  xx=np.arange(a,b)*dt; yy=y[a:b];xc=xx-xx.mean(); den=xc@xc
  if den>0:out[i]=max(0,(xc@(yy-yy.mean()))/den)
 return pd.Series(out).interpolate(limit_direction='both').to_numpy()
def smooth_deriv(y,dt,sec=5):
 w=max(5,int(round(sec/dt)));w+=1-(w%2);w=min(w,len(y) if len(y)%2 else len(y)-1)
 ys=savgol_filter(y,w,2,mode='interp') if w>=5 else y
 return np.gradient(ys,dt)
def corr_profile(x,y,dt,maxlag=60):
 n=min(len(x),len(y));x=x[:n];y=y[:n];K=min(int(round(maxlag/dt)),n-5)
 cr=correlate(y,x,mode='full',method='fft')[n-1:n+K]
 px=np.r_[0,np.cumsum(x)];py=np.r_[0,np.cumsum(y)];px2=np.r_[0,np.cumsum(x*x)];py2=np.r_[0,np.cumsum(y*y)]
 k=np.arange(K+1);m=n-k;sx=px[n-k];sy=py[n]-py[k];sx2=px2[n-k];sy2=py2[n]-py[k]
 cov=cr-sx*sy/m;den=np.sqrt(np.maximum(0,sx2-sx*sx/m)*np.maximum(0,sy2-sy*sy/m))
 r=np.divide(cov,den,out=np.full_like(cov,np.nan),where=den>1e-12)
 return k*dt,np.clip(r,-1,1)

for gas in GASES:
 fig,axs=plt.subplots(2,5,figsize=(18,7),sharex=True,constrained_layout=True)
 for ax,h in zip(axs.flat,HEATS):
  d=pd.read_csv(ROOT/'inputs'/f'curves_{h}.csv');dt=.5
  q=roll_slope(d.oxygen_cumulative.to_numpy(float),dt);x=smooth_deriv(q,dt);y=smooth_deriv(d[f'{gas}_flow'].to_numpy(float),dt)
  trim=int(round(20/dt));x=x[trim:-trim];y=y[trim:-trim];lag,r=corr_profile(x,y,dt)
  ax.plot(lag,r,lw=1.2);ax.axhline(0,lw=.6);ax.set_title(h,fontsize=9);ax.grid(True,alpha=.2)
  ax.set_xlabel('Lag (s)');ax.set_ylabel('Pearson r')
 fig.suptitle(f'Oxygen-flow change → {gas} change: H0 full-series lag profiles',fontsize=14)
 fig.savefig(OUT/f'figures/H0_full_lag_profiles_{gas}.png',dpi=200,bbox_inches='tight');plt.close(fig)

for scenario in ['H0','H1']:
 fig,axs=plt.subplots(1,2,figsize=(14,5.5),constrained_layout=True)
 for ax,gas in zip(axs,GASES):
  sub=windows[(windows.scenario==scenario)&(windows.gas==gas)]
  for h,g in sub.groupby('heat'): ax.scatter(g.flow_median_m3s,g.best_lag_s,s=25,alpha=.7,label=h)
  ax.axvspan(windows[windows.heat!='6C02057'].flow_median_m3s.quantile(.01),windows[windows.heat!='6C02057'].flow_median_m3s.quantile(.99),alpha=.08)
  ax.set_xlabel('Median oxygen flow (m³/s)');ax.set_ylabel('Local lag (s)');ax.set_title(gas);ax.grid(True,alpha=.25);ax.legend(fontsize=7,ncol=2,frameon=False)
 fig.suptitle(f'Flow vs local lag — {scenario}',fontsize=14);fig.savefig(OUT/f'figures/{scenario}_flow_vs_local_lag.png',dpi=200,bbox_inches='tight');plt.close(fig)

for scenario in ['H0','H1']:
 fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
 for ax,gas in zip(axs,GASES):
  z=bins[(bins.scenario==scenario)&(bins.gas==gas)&(bins.window_set=='informative')].set_index('flow_bin').reindex(['low','middle','high'])
  x=np.arange(3);med=z.lag_median_s.to_numpy(float);lo=med-z.lag_q25_s.to_numpy(float);hi=z.lag_q75_s.to_numpy(float)-med
  ax.errorbar(x,med,yerr=[np.maximum(lo,0),np.maximum(hi,0)],fmt='o-',capsize=4)
  ax.set_xticks(x,labels=['Low','Middle','High']);ax.set_ylabel('Lag (s)');ax.set_title(gas);ax.grid(True,alpha=.25)
 fig.suptitle(f'Flow-bin lag medians and IQR — {scenario}',fontsize=14);fig.savefig(OUT/f'figures/{scenario}_flow_bin_lag.png',dpi=200,bbox_inches='tight');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
for ax,gas in zip(axs,GASES):
 z=full[(full.heat=='6C02057')&(full.gas==gas)].set_index('scenario').reindex(['H0','H1'])
 ax.bar(['H0 0.5s','H1 1.0s'],z.best_abs_lag_s.to_numpy());ax.set_ylabel('Best |r| lag (s)');ax.set_title(gas);ax.grid(axis='y',alpha=.25)
fig.suptitle('Heat 6C02057: time-scale sensitivity',fontsize=14);fig.savefig(OUT/'figures/6C02057_H0_H1_lag_sensitivity.png',dpi=200,bbox_inches='tight');plt.close(fig)

if len(spec):
 fig,axs=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
 for ax,gas in zip(axs,GASES):
  z=spec[(spec.scenario=='H0')&(spec.gas==gas)]
  for sw,g in z.groupby('smooth_s'):ax.scatter(g.best_abs_lag_s,g.best_abs_r.abs(),label=f'{sw:g}s',alpha=.7)
  ax.set_xlabel('Best absolute lag (s)');ax.set_ylabel('Peak |r|');ax.set_title(gas);ax.grid(True,alpha=.25);ax.legend(frameon=False)
 fig.suptitle('Supplementary optical brightness → gas lag diagnostic (H0)',fontsize=14);fig.savefig(OUT/'figures/H0_spectral_to_gas_lag_diagnostic.png',dpi=200,bbox_inches='tight');plt.close(fig)

# Summary and report.
summary={}
for sc in ['H0','H1']:
 summary[sc]={}
 for gas in GASES:
  z=full[(full.scenario==sc)&(full.gas==gas)]
  w=windows[(windows.scenario==sc)&(windows.gas==gas)];wi=w[w.informative]
  summary[sc][gas]={'median_peak_lag_s':float(z.best_abs_lag_s.median()),
   'median_peak_abs_r':float(z.best_abs_r.abs().median()),
   'BH_q_le_0_10_count':int((z.BH_q<=.10).sum()),
   'median_bootstrap_ci_width_s':float((z.bootstrap_lag_ci_high_s-z.bootstrap_lag_ci_low_s).median()),
   'informative_windows':len(wi),'informative_heats':wi.heat.nunique(),
   'informative_lag_median_s':float(wi.best_lag_s.median()) if len(wi) else None,
   'informative_flow_median_m3s':float(wi.flow_median_m3s.median()) if len(wi) else None}
summary['flow_slopes']=slopes.replace({np.nan:None}).to_dict(orient='records')
summary['loho_model_mae']=models.groupby(['scenario','gas','model']).test_MAE_s.mean().reset_index().to_dict(orient='records') if len(models) else []
summary['spectral_H0']={}
if len(spec):
 for gas in GASES:
  z=spec[(spec.scenario=='H0')&(spec.gas==gas)]
  summary['spectral_H0'][gas]={'median_peak_abs_r':float(z.best_abs_r.abs().median()),'median_lag_s':float(z.best_abs_lag_s.median()),'boundary_peaks':int(z.boundary.sum()),'n':len(z)}
(OUT/'tables/summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

# README report with concise numerical findings.
report=['# 第二问延迟识别专项实验（Q2-delay-v1.0）','',
'**数据口径：** 十炉曲线；H0 全部 0.5 s/行；H1 仅 6C02057 按 1.0 s/行。H0 为官方主情景，H1 为敏感性情景。', '',
'## 实验问题','是否存在可重复的供氧变化→烟气变化有效时差，以及该时差是否随供氧流量增加而缩短。', '',
'## 方法摘要','- 由累计耗氧量在 10 s 居中窗口内局部线性回归得到非负供氧流量估计；以氧流和气体流量各自 5 s 平滑后的一阶变化率计算互相关。正滞后表示气体响应较晚。全炉搜索 0–60 s。',
'- 120 s 窗口、60 s 步长、局部搜索 0–40 s；信息性窗口定义为氧流标准差≥0.15 m³/s、峰值 |r|≥0.25、峰不在搜索边界。该阈值仅为探索性筛选。',
'- 全炉最大相关峰用 120 s 配对移动块 bootstrap 估计峰值滞后的 95%区间；循环移位作为保持单序列自相关的近似负对照，BH 校正跨 20 个炉次×气体检验。',
'- 流量—延迟斜率使用炉内中心化，避免炉间固定差异伪装成流量效应；按炉次聚类 bootstrap。H0 的通用斜率不纳入高流量第五炉外推点；H1 则第五炉回到九炉共同流量支持范围附近。',
'- 补充计算 32 波段对数亮度代理到 CO₂/CO 的延迟。它不是第三问最终冻结的光谱状态。','']
report += ['## 结果速览','']
for sc in ['H0','H1']:
 for gas in GASES:
  s=summary[sc][gas]
  report.append(f"- {sc} / {gas}：全炉峰值时差中位数 {s['median_peak_lag_s']:.1f} s；峰值 |r| 中位数 {s['median_peak_abs_r']:.3f}；BH q≤0.10 的炉次 {s['BH_q_le_0_10_count']}/10；bootstrap 区间宽度中位数 {s['median_bootstrap_ci_width_s']:.1f} s；信息性局部窗口 {s['informative_windows']} 个、覆盖 {s['informative_heats']} 炉。")
report.append('')
report.append('### 流量—延迟斜率')
for r in slope_rows:
 report.append(f"- {r['scenario']} / {r['gas']}：炉内中心化斜率 {r['within_heat_slope_s_per_m3s']:.3f} s/(m³/s)，炉次聚类 bootstrap 95%区间 [{r['cluster_boot_ci_low']:.3f}, {r['cluster_boot_ci_high']:.3f}]，负斜率重抽比例 {r['bootstrap_share_beta_negative']:.3f}，LOHO 负斜率比例 {r['loo_negative_share']:.3f}。")
report.append('')
report.append('### 固定 vs 流量条件模型（留一炉）')
if len(models):
 mm=models.groupby(['scenario','gas','model']).test_MAE_s.mean()
 for idx,v in mm.items():report.append(f'- {idx[0]} / {idx[1]} / {idx[2]}：MAE {v:.2f} s。')
report += ['', '## 解释边界',
'- 相关峰不是物理运输时间真值；炉内反应、气体混合、烟道输运、仪器响应和时钟偏差无法拆分。',
'- 氧枪瞬时流量是从累计耗氧反推，没有独立流量传感器标签。',
'- 10 炉样本少；同炉窗口重叠且自相关强，窗口数量不能当作独立炉次数。',
'- 循环移位是近似负对照，非平稳序列下不是严格随机化检验；bootstrap 区间也受弱激励与峰值选择影响。',
'- 不能因题面写“流量越大、延迟越小”就强制非递增模型；只有负斜率跨炉稳定、区间排除零、LOHO 改善时才可采用。',
'- H0/H1 的第五炉差异应并列披露；当前数据不能单独识别其绝对采样时间尺度。', '',
'## 关键文件','- `tables/full_series_lag_diagnostics.csv`：全炉延迟峰、峰相关、bootstrap 区间、循环移位 p 与 BH q。',
'- `tables/local_window_lag_diagnostics.csv`：逐窗口流量和延迟诊断。',
'- `tables/flow_bin_lag_summary.csv`：流量分箱。',
'- `tables/within_heat_flow_delay_slope.csv`：流量—延迟斜率及炉级 bootstrap 区间。',
'- `tables/flow_delay_model_loho.csv`：固定、自由线性、非递增等模型的留一炉误差。',
'- `tables/spectral_brightness_lag_diagnostics.csv`：辅助光谱亮度—烟气时差。',
'- `tables/gas_signal_quality.csv`：气体序列质量摘要。',
'- `code/run_q2_delay_experiment.py`：主实验脚本（首次运行的原始计算脚本）。',
'- `code/postprocess_q2_delay.py`：结果汇总与图表脚本。', '',
'## 当前判断','若流量—延迟斜率的炉级区间跨零、符号不稳定，或非递增模型的留一炉误差未优于固定基准，则当前数据不足以支持普遍的“流量越大、延迟越小”定量规律。此时应报告有效时差的探索性范围和弱识别结论，不虚构精确参数。']
(OUT/'README_Q2_Delay_Experiment.md').write_text('\n'.join(report),encoding='utf-8')

# Save this postprocess script.
shutil.copy2(Path(__file__),OUT/'code/postprocess_q2_delay.py')

zip_path=ROOT.parent/'Q2_Delay_Identification_Experiment_rerun.zip'
if zip_path.exists():zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():zf.write(p,arcname=str(p.relative_to(OUT)))

print('ZIP',zip_path,zip_path.stat().st_size)
print('\nFULL')
print(full[['scenario','heat','gas','best_abs_lag_s','best_abs_r','bootstrap_lag_ci_low_s','bootstrap_lag_ci_high_s','circular_shift_p','BH_q','boundary']].round(3).to_string(index=False))
print('\nSLOPES')
print(slopes.round(3).to_string(index=False))
print('\nBINS informative')
print(bins[bins.window_set=='informative'].round(2).to_string(index=False))
print('\nLOHO')
print(models.groupby(['scenario','gas','model']).test_MAE_s.mean().round(2).to_string() if len(models) else 'no models')
print('\nSPECTRAL')
print(spec[spec.scenario=='H0'].groupby(['gas','smooth_s']).agg(median_abs_r=('best_abs_r',lambda x:np.median(np.abs(x))),median_lag=('best_abs_lag_s','median'),boundary=('boundary','sum')).round(3).to_string() if len(spec) else 'none')
