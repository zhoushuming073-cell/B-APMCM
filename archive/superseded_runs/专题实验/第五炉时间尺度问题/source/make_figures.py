from pathlib import Path
import argparse, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import welch
sys.path.insert(0,str(Path(__file__).resolve().parent)); import analysis_helpers as base
parser=argparse.ArgumentParser(); parser.add_argument('--cache-dir',default='cache'); parser.add_argument('--out',default='results'); args=parser.parse_args()
CACHE=Path(args.cache_dir).resolve(); OUT=Path(args.out).resolve(); FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
D={}
for f in CACHE.glob('heats_*.npz'):
 z=np.load(f); h=f.stem.replace('heats_',''); D[h]={'seq':z['seq'],'X':z['X'],'proc':z['proc'],'wl':z['wl']}
TARGET='6C02057';
# Features
F={}
for h,d in D.items(): F[h]=base.spectral_features(d['X'],d['wl'])
P= D[TARGET]['proc']; O=P[:,0]
flow=pd.read_csv(OUT/'tables'/'oxygen_flow_summary.csv')
lag=pd.read_csv(OUT/'tables'/'lag_estimates.csv')
lagmod=pd.read_csv(OUT/'tables'/'flow_conditioned_lag_models.csv')
lagcoef=pd.read_csv(OUT/'tables'/'lag_model_coefficients_and_summary.csv')
macro=pd.read_csv(OUT/'tables'/'distributed_lag_LOHO_fair_macro.csv')
gasint=pd.read_csv(OUT/'tables'/'gas_integrals_and_balance.csv')
quant=pd.read_csv(OUT/'tables'/'counter_quantization_diagnostics.csv')
cv=pd.read_csv(OUT/'tables'/'fifth_oxygen_spline_cv.csv')
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'figure.dpi':160})
# Fig 01: all row differences including zeros
fig,axs=plt.subplots(2,5,figsize=(15,6),constrained_layout=True)
for ax,(h,d) in zip(axs.flat,D.items()):
 delta=np.diff(d['proc'][:,0]); bins=np.r_[-.5,.5,np.arange(1,int(max(np.max(np.diff(x['proc'][:,0])) for x in D.values()))+2,1)]
 ax.hist(delta,bins=bins,density=True,alpha=.8)
 ax.set_title(f'{h}\nzero={np.mean(delta==0):.1%}'); ax.set_xlabel('Δ cumulative O₂ (m³)'); ax.set_ylabel('Density')
fig.suptitle('All one-row cumulative-oxygen differences (zeros retained)'); fig.savefig(FIG/'fig01_oxygen_delta_distributions.png'); plt.close(fig)
# Fig 02: update gap distribution
fig,axs=plt.subplots(2,5,figsize=(15,6),constrained_layout=True)
for ax,(h,d) in zip(axs.flat,D.items()):
 delta=np.diff(d['proc'][:,0]); ix=np.flatnonzero(delta>0)+1; g=np.diff(np.r_[0,ix]); bins=np.arange(.5,max(8.5,g.max()+1.5),1)
 ax.hist(g,bins=bins,density=True,alpha=.8); ax.axvline(np.median(g),ls='--',lw=1)
 ax.set_title(f'{h}\nmedian={np.median(g):.0f} rows'); ax.set_xlabel('Rows between updates')
fig.suptitle('Counter update gaps measured in rows (not sampling intervals)'); fig.savefig(FIG/'fig02_update_gap_distributions.png'); plt.close(fig)
# Fig 03 cumulative O2
fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
for sc,dt in [('H0',.5),('H1',1.)]: ax.plot(np.arange(len(O))*dt,O,lw=1.6,label=f'{sc}: Δt={dt:g} s')
ax.set_xlabel('Elapsed time under hypothesis (s)'); ax.set_ylabel('Cumulative oxygen (m³)'); ax.legend(); ax.grid(alpha=.25); ax.set_title('6C02057 cumulative oxygen: identical observations, different clocks'); fig.savefig(FIG/'fig03_fifth_cumulative_oxygen_H0_H1.png'); plt.close(fig)
# Fig 04 flow curve
fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
for sc,dt in [('H0',.5),('H1',1.)]:
 q=base.flow_slope(O,dt,10.,False); ax.plot(np.arange(len(O))*dt,q,label=f'{sc}: 10-s centered slope')
ax.set_xlabel('Elapsed time under hypothesis (s)'); ax.set_ylabel('Estimated oxygen flow (m³/s)'); ax.legend(); ax.grid(alpha=.25); ax.set_title('6C02057 oxygen flow under H0/H1'); fig.savefig(FIG/'fig04_fifth_flow_H0_H1.png'); plt.close(fig)
# Fig 05 spectral six-panel dynamics
fig,axs=plt.subplots(3,2,figsize=(13,9),constrained_layout=True)
items=[('spectral_integral','Integrated spectral intensity'),('spectral_p95','Spectral 95th percentile'),('spectral_centroid','Spectral centroid (nm)'),('peak1_height','Dominant peak height'),('peak1_width_nm','Dominant peak width (nm)'),('peak1_area','Dominant peak area')]
for ax,(key,label) in zip(axs.flat,items):
 for sc,dt in [('H0',.5),('H1',1.)]: ax.plot(np.arange(len(O))*dt,F[TARGET][key],lw=1,label=sc)
 ax.set_title(label); ax.set_xlabel('Elapsed time (s)'); ax.grid(alpha=.2); ax.legend()
fig.suptitle('6C02057 raw spectral feature dynamics under H0/H1'); fig.savefig(FIG/'fig05_fifth_spectral_dynamics.png'); plt.close(fig)
# Fig 06 gas curves
fig,axs=plt.subplots(2,1,figsize=(11,7),constrained_layout=True)
for sc,dt in [('H0',.5),('H1',1.)]:
 t=np.arange(len(O))*dt; axs[0].plot(t,P[:,1],label=sc); axs[1].plot(t,P[:,2],label=sc)
axs[0].set_ylabel('CO₂ volume flow (m³/s)'); axs[1].set_ylabel('CO volume flow (m³/s)'); axs[1].set_xlabel('Elapsed time (s)')
for ax in axs: ax.grid(alpha=.2); ax.legend()
fig.suptitle('6C02057 gas-flow signals under H0/H1'); fig.savefig(FIG/'fig06_fifth_gas_dynamics.png'); plt.close(fig)
# Fig 07 spectral ACF
fig,axs=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
for j,key in enumerate(['spectral_integral','peak1_height']):
 for h in D:
  if h==TARGET or h.startswith('6B'): continue
  a=base.acf(F[h][key],150); axs[j].plot(np.arange(len(a))*.5,a,color='0.75',alpha=.35,lw=.8)
 for sc,dt in [('H0',.5),('H1',1.)]:
  a=base.acf(F[TARGET][key],150); axs[j].plot(np.arange(len(a))*dt,a,lw=2,label=sc)
 axs[j].set_title(key); axs[j].set_xlabel('Lag (s)'); axs[j].set_ylabel('ACF'); axs[j].grid(alpha=.2); axs[j].legend()
fig.suptitle('Spectral ACF: other 6C heats at official 0.5 s'); fig.savefig(FIG/'fig07_spectral_acf.png'); plt.close(fig)
# Fig 08 gas ACF
fig,axs=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
for j,key in enumerate(['CO2','CO']):
 for h,d in D.items():
  if h==TARGET or h=='6B03205': continue
  a=base.acf(d['proc'][:,j+1],150); axs[j].plot(np.arange(len(a))*.5,a,color='0.75',alpha=.35,lw=.8)
 for sc,dt in [('H0',.5),('H1',1.)]:
  a=base.acf(P[:,j+1],150); axs[j].plot(np.arange(len(a))*dt,a,lw=2,label=sc)
 axs[j].set_title(key); axs[j].set_xlabel('Lag (s)'); axs[j].set_ylabel('ACF'); axs[j].grid(alpha=.2); axs[j].legend()
fig.suptitle('Gas-signal ACF: other heats at official 0.5 s'); fig.savefig(FIG/'fig08_gas_acf.png'); plt.close(fig)
# Fig 09 fair model MAE
fig,axs=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
for j,g in enumerate(['CO2','CO']):
 p=macro[macro.gas==g]; names=['M0_persistence','M1_current_flow','M2_distributed_lag','M3_history_windows']; x=np.arange(len(names)); w=.36
 for k,sc in enumerate(['H0','H1']):
  ss=p[p.scenario==sc].set_index('model').reindex(names); axs[j].bar(x+(k-.5)*w,ss.MAE_macro,w,label=sc)
 axs[j].set_xticks(x,names,rotation=20,ha='right'); axs[j].set_ylabel('LOHO macro MAE (m³/s)'); axs[j].set_title(g); axs[j].grid(axis='y',alpha=.2); axs[j].legend()
fig.suptitle('Matched-row LOHO prediction MAE under H0/H1'); fig.savefig(FIG/'fig09_model_MAE_H0_H1.png'); plt.close(fig)
# Fig 10 lag estimates
fig,axs=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
heats=[h for h in D if h!='6B03205']
for j,g in enumerate(['CO2','CO']):
 for sc,clock,color,marker in [('H0','H0','C0','o'),('H1','H1','C1','s')]:
  others=lag[(lag.heat!=TARGET)&(lag.heat!='6B03205')&(lag.clock=='official')&(lag.gas==g)]
  target=lag[(lag.heat==TARGET)&(lag.clock==clock)&(lag.gas==g)]
  ss=pd.concat([others,target]); axs[j].scatter([heats.index(h) for h in ss.heat],ss.lag_s,label=sc,color=color,marker=marker)
 axs[j].set_xticks(range(len(heats)),heats,rotation=45,ha='right'); axs[j].set_ylabel('Estimated lag (s)'); axs[j].set_title(g); axs[j].grid(alpha=.2); axs[j].legend()
fig.suptitle('Exploratory flow-to-gas lag estimates; cross-correlation only'); fig.savefig(FIG/'fig10_lag_estimates.png'); plt.close(fig)
# Fig 11 flow-lag scatter
fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
for j,g in enumerate(['CO2','CO']):
 for sc,clock,color in [('H0','H0','C0'),('H1','H1','C1')]:
  others=lag[(lag.heat!=TARGET)&(lag.heat!='6B03205')&(lag.clock=='official')&(lag.gas==g)]
  target=lag[(lag.heat==TARGET)&(lag.clock==clock)&(lag.gas==g)]
  ss=pd.concat([others,target]); axs[j].scatter(ss.flow_median_m3_s,ss.lag_s,label=sc,color=color)
  for _,r in ss.iterrows(): axs[j].annotate(r.heat[-5:],(r.flow_median_m3_s,r.lag_s),fontsize=6,xytext=(2,2),textcoords='offset points')
 axs[j].set_xlabel('Median oxygen flow (m³/s)'); axs[j].set_ylabel('Estimated lag (s)'); axs[j].set_title(g); axs[j].grid(alpha=.2); axs[j].legend()
fig.suptitle('Flow-conditioned delay scatter under H0/H1'); fig.savefig(FIG/'fig11_flow_lag_sensitivity.png'); plt.close(fig)
# Fig 12 gas integrated volumes
fig,axs=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
for j,g in enumerate(['CO2_plus_CO','CO2','CO']):
 if j>1: break
 sub=gasint[gasint.gas==g]
 other=sub[sub.heat!=TARGET].set_index('heat').integral_m3
 vals=list(other.values)
 target=sub[sub.heat==TARGET].set_index('clock').integral_m3
 axs[j].boxplot(vals,positions=[1],widths=.5); axs[j].scatter([1],[target['H0']],label='5th H0',marker='o'); axs[j].scatter([1],[target['H1']],label='5th H1',marker='s')
 axs[j].set_xticks([1],['Other 9 heats']); axs[j].set_ylabel('Integrated volume (m³)'); axs[j].set_title(g); axs[j].legend(); axs[j].grid(axis='y',alpha=.2)
fig.suptitle('Gas-flow integrals: H1 doubles the fifth-heat totals'); fig.savefig(FIG/'fig12_gas_integral_consistency.png'); plt.close(fig)
# Fig 13 quantization
fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
heats=list(quant.heat); x=np.arange(len(heats));
axs[0].bar(x,quant.positive_delta_mean,label='Mean ΔO'); axs[0].scatter(x,quant.positive_delta_median,color='black',label='Median ΔO',zorder=3)
axs[0].set_xticks(x,heats,rotation=45,ha='right'); axs[0].set_ylabel('m³ per positive update'); axs[0].set_title('Mean vs median positive increments'); axs[0].legend(); axs[0].grid(axis='y',alpha=.2)
axs[1].bar(x,quant.share_increment_2plus_quanta*100); axs[1].set_xticks(x,heats,rotation=45,ha='right'); axs[1].set_ylabel('Share of positive updates (%)'); axs[1].set_title('Positive increments ≥2 × 4.74 m³'); axs[1].grid(axis='y',alpha=.2)
fig.suptitle('The fifth-heat upper tail is hidden by the median'); fig.savefig(FIG/'fig13_counter_increment_quantization.png'); plt.close(fig)
# Fig 14 update gap vs increment target
z=D[TARGET]; delta=np.diff(z['proc'][:,0]); ix=np.flatnonzero(delta>0)+1; gaps=np.diff(np.r_[0,ix]); p=delta[delta>0]
fig,ax=plt.subplots(figsize=(8,5),constrained_layout=True)
ax.scatter(gaps,p,alpha=.45,s=18); ax.set_xlabel('Rows between counter updates'); ax.set_ylabel('Positive ΔO (m³)'); ax.set_title('6C02057: update gap and increment are weakly associated'); ax.grid(alpha=.2)
fig.savefig(FIG/'fig14_fifth_update_gap_vs_increment.png'); plt.close(fig)
# Fig 15 PSD of spectral integral
fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
for h in D:
 if h==TARGET or h.startswith('6B'): continue
 x=F[h]['spectral_integral']; f,p=welch(x-np.mean(x),fs=2.,nperseg=256); ax.plot(f,p,alpha=.35,color='0.65',lw=.8)
for sc,dt,color in [('H0',.5,'C0'),('H1',1.,'C1')]:
 x=F[TARGET]['spectral_integral']; f,p=welch(x-np.mean(x),fs=1/dt,nperseg=256); ax.plot(f,p,lw=2,color=color,label=sc)
ax.set_xlim(0,0.5); ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel('PSD'); ax.set_title('Spectral-integral PSD: other 6C heats at 0.5 s'); ax.grid(alpha=.2); ax.legend()
fig.savefig(FIG/'fig15_spectral_PSD.png'); plt.close(fig)
# Fig 16 adjacent spectrum similarity
fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
vals=[]; labels=[]
for h in D:
 vals.append(F[h]['adjacent_cosine'][np.isfinite(F[h]['adjacent_cosine'])]); labels.append(h)
ax.boxplot(vals,labels=labels,showfliers=False); ax.set_xticklabels(labels,rotation=45,ha='right'); ax.set_ylabel('Cosine similarity'); ax.set_title('Adjacent-row spectral similarity (index-based; invariant to clock scaling)'); ax.grid(axis='y',alpha=.2)
fig.savefig(FIG/'fig16_adjacent_spectral_similarity.png'); plt.close(fig)
# Fig 17 flow median comparison
fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
sub=flow[flow.method=='centered_slope_10s']; heats=[h for h in D]; x=np.arange(len(heats)); w=.36
for k,clock in enumerate(['official','H0','H1']):
 ss=sub[sub.clock==clock].set_index('heat').reindex(heats)
 if ss['p50'].notna().any(): ax.bar(x+(k-1)*w,ss.p50,w,label=clock)
ax.set_xticks(x,heats,rotation=45,ha='right'); ax.set_ylabel('10-s centered flow median (m³/s)'); ax.set_title('Cross-heat oxygen-flow medians under H0/H1'); ax.legend(); ax.grid(axis='y',alpha=.2)
fig.savefig(FIG/'fig17_flow_median_by_heat.png'); plt.close(fig)
# Fig 18 flow spline CV
fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
ax.bar(cv.clock,cv.integral_spline_outer_MAE); ax.set_ylabel('Outer-block MAE (m³)'); ax.set_title('6C02057 integral-spline cumulative-O₂ fit'); ax.grid(axis='y',alpha=.2)
fig.savefig(FIG/'fig18_oxygen_spline_CV.png'); plt.close(fig)
print('wrote',len(list(FIG.glob('*.png'))),'figures to',FIG)
