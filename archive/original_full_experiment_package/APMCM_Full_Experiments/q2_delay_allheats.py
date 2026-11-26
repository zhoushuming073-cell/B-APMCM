import os,glob
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
R='/mnt/data/APMCM_Full_Experiments';files=sorted(glob.glob('/mnt/data/heats_*.npz'))
D={}
for f in files:
 z=np.load(f);h=os.path.basename(f)[6:-4];D[h]={k:z[k] for k in z.files}

def flow_slope(O,dt=.5,win=10,causal=False):
 n=len(O);t=np.arange(n)*dt;k=int(round(win/dt));q=np.full(n,np.nan)
 for i in range(n):
  a=max(0,i-k+1) if causal else max(0,i-k//2);b=i+1 if causal else min(n,i+k//2+1)
  if b-a<3:continue
  x=t[a:b];y=O[a:b];xc=x-x.mean();q[i]=max(0,xc@(y-y.mean())/(xc@xc))
 return q

def corr_lags(x,y,maxlag=80):
 # positive lag means gas change follows upstream flow change
 out=[]
 for lag in range(maxlag+1):
  a=x[:len(x)-lag] if lag else x
  b=y[lag:] if lag else y
  ok=np.isfinite(a)&np.isfinite(b)
  if ok.sum()<50 or np.std(a[ok])<1e-8 or np.std(b[ok])<1e-8:r=np.nan
  else:r=float(np.corrcoef(a[ok],b[ok])[0,1])
  out.append(r)
 return np.array(out)

lag_rows=[];window_rows=[];corr_rows=[];curves={}
for h,z in D.items():
 n=len(z['seq']);dt=.5;t=np.arange(n)*dt;O=z['proc'][:,0].astype(float);G=z['proc'][:,1:3].astype(float);qcen=flow_slope(O,win=10,causal=False);qtrail=flow_slope(O,win=10,causal=True)
 # retrospective cross-correlation on smoothed first differences
 qsm=uniform_filter1d(np.nan_to_num(qcen,nan=np.nanmedian(qcen)),size=10,mode='nearest')
 dq=np.diff(qsm,prepend=qsm[0])
 for j,gn in enumerate(['CO2','CO']):
  y=G[:,j];ysm=uniform_filter1d(y,size=10,mode='nearest');dy=np.diff(ysm,prepend=ysm[0]);cc=corr_lags(dq,dy,80);valid=np.isfinite(cc)
  if valid.any():
   best=int(np.nanargmax(np.abs(cc)));bestpos=int(np.nanargmax(cc));bestneg=int(np.nanargmin(cc))
   # Circular-shift negative control of derivative series; exclude small shifts to preserve no near-zero alignment.
   rng=np.random.default_rng(20261126+sum(map(ord,h+gn)));obs=np.nanmax(np.abs(cc));null=[]
   for s in rng.integers(80,max(81,n-80),size=199):
    null.append(np.nanmax(np.abs(corr_lags(dq,np.roll(dy,int(s)),80))))
   p=(1+sum(v>=obs for v in null))/200
   lag_rows.append({'heat':h,'gas':gn,'best_abs_lag_s':best*dt,'best_abs_corr':float(cc[best]),'best_positive_lag_s':bestpos*dt,'best_positive_corr':float(cc[bestpos]),'best_negative_lag_s':bestneg*dt,'best_negative_corr':float(cc[bestneg]),'shift_p':p,'q_sd':float(np.nanstd(qcen)),'q_p10':float(np.nanquantile(qcen,.1)),'q_p50':float(np.nanquantile(qcen,.5)),'q_p90':float(np.nanquantile(qcen,.9)),'d_q_sd':float(np.std(dq)),'status':'exploratory_cross_correlation'})
  # 60s nonoverlap windows, report only; no forced positive delay
  for st in range(0,int(t[-1])-60,60):
   a=int(st/dt);b=min(n,int((st+60)/dt));x=dq[a:b];yy=dy[a:b]
   if len(x)<80:continue
   ccw=corr_lags(x,yy,80);k=int(np.nanargmax(abs(ccw))) if np.any(np.isfinite(ccw)) else -1
   window_rows.append({'heat':h,'gas':gn,'window_start_s':st,'flow_median':float(np.nanmedian(qcen[a:b])),'flow_sd':float(np.nanstd(qcen[a:b])),'best_lag_s':k*dt if k>=0 else np.nan,'best_corr':float(ccw[k]) if k>=0 else np.nan,'boundary':bool(k in [0,80]) if k>=0 else True,'informative':bool(k>=0 and np.nanstd(qcen[a:b])>=.15 and not (k in [0,80]) and abs(ccw[k])>=.3)})
 # associations and save curve
 cO=np.r_[0,np.cumsum((G[:-1,0]+G[1:,0])*.25)];cCO=np.r_[0,np.cumsum((G[:-1,1]+G[1:,1])*.25)]
 curves[h]=pd.DataFrame({'time_s':t,'oxygen_cumulative':O,'oxygen_flow_centered10':qcen,'oxygen_flow_trailing10':qtrail,'CO2_flow':G[:,0],'CO_flow':G[:,1],'CO2_cumulative':cO,'CO_cumulative':cCO})
 print('delay done',h,flush=True)
# LOHO forward prediction at target time t, issue time t-5s; 1Hz downsample.
# Feature sets: persistence + process progress; then add q current/history/lagged samples.
heat_arrays={}
for h,z in D.items():
 n=len(z['seq']);idx=np.arange(0,n,2);t=idx*.5;O=z['proc'][:,0][idx];G=z['proc'][:,1:3][idx];q=flow_slope(z['proc'][:,0].astype(float),win=10,causal=True)[idx]
 q=np.nan_to_num(q,nan=np.nanmedian(q))
 # past means over intervals preceding issue time (1Hz index)
 def pastmean(x,k,a,b):
  out=np.zeros(len(x));
  for i in range(len(x)):
   en=max(0,i-k);st=max(0,i-b);en2=max(0,i-a)
   out[i]=np.mean(x[st:en2]) if en2>st else x[max(0,i-k)]
  return out
 # Target index i; issue at i-5s; all q and gas features use <= issue.
 n1=len(idx);issue=np.arange(n1)-5;target=np.arange(n1)
 valid=(issue>=40)&(target>=45)
 def atlag(x,lag):return x[np.maximum(0,issue-lag)]
 prog=np.c_[t/max(t[-1],1),O/max(O[-1],1)]
 pers=np.column_stack([G[np.maximum(0,issue),0],G[np.maximum(0,issue),1]])
 qlag=np.column_stack([atlag(q,l) for l in [0,5,10,20,30,40]])
 qhist=np.column_stack([np.array([np.mean(q[max(0,k-a-b):max(0,k-a)]) if max(0,k-a)>max(0,k-a-b) else q[max(0,k-a)] for k in issue]) for a,b in [(0,5),(5,10),(15,15),(30,10),(40,20)]])
 # Family dummy controls optical/data source group; no calibration family indicator as flow data are not optical.
 heat_arrays[h]={'valid':valid,'prog':prog[np.maximum(0,issue)],'pers':pers,'qcurrent':qlag[:,:1],'qhist':qhist,'qlag':qlag,'G':G[target],'t':t[target]}
# outer leave-one-heat-out, inner LOHO alpha selection
pred_rows=[];models={'M0_persistence':lambda d: d['pers'],'M1_progress_persistence':lambda d: np.c_[d['pers'],d['prog']],'M2_current_flow':lambda d: np.c_[d['pers'],d['prog'],d['qcurrent']],'M3_history_flow':lambda d: np.c_[d['pers'],d['prog'],d['qhist']],'M4_lag_grid':lambda d: np.c_[d['pers'],d['prog'],d['qlag']]}
for testh in heat_arrays:
 trainhs=[h for h in heat_arrays if h!=testh]
 for j,gn in enumerate(['CO2','CO']):
  # same rows for every model
  train_idx={h:np.where(heat_arrays[h]['valid'])[0] for h in trainhs};test_idx=np.where(heat_arrays[testh]['valid'])[0]
  ytest=heat_arrays[testh]['G'][test_idx,j]
  for mn,fn in models.items():
   # tune alpha with inner leave-one-heat-out macro MAE
   alphas=[.1,1,10,100];inner=[]
   for alpha in alphas:
    vals=[]
    for vh in trainhs:
     itrain=[h for h in trainhs if h!=vh]
     if not itrain:continue
     Xtr=np.vstack([fn(heat_arrays[h])[train_idx[h]] for h in itrain]);ytr=np.concatenate([heat_arrays[h]['G'][train_idx[h],j] for h in itrain])
     Xv=fn(heat_arrays[vh])[train_idx[vh]];yv=heat_arrays[vh]['G'][train_idx[vh],j]
     model=make_pipeline(StandardScaler(),Ridge(alpha=alpha));model.fit(Xtr,ytr);yp=np.maximum(0,model.predict(Xv));vals.append(np.mean(abs(yp-yv)))
    inner.append(np.mean(vals))
   alpha=alphas[int(np.argmin(inner))]
   Xtr=np.vstack([fn(heat_arrays[h])[train_idx[h]] for h in trainhs]);ytr=np.concatenate([heat_arrays[h]['G'][train_idx[h],j] for h in trainhs]);Xte=fn(heat_arrays[testh])[test_idx]
   model=make_pipeline(StandardScaler(),Ridge(alpha=alpha));model.fit(Xtr,ytr);yp=np.maximum(0,model.predict(Xte))
   pred_rows.append({'test_heat':testh,'gas':gn,'model':mn,'alpha':alpha,'inner_mae':float(min(inner)),'test_MAE':float(np.mean(abs(yp-ytest))),'test_RMSE':float(np.sqrt(np.mean((yp-ytest)**2))),'n_test':len(ytest)})
# Save outputs
pd.DataFrame(lag_rows).to_csv(R+'/q2_lag_by_heat.csv',index=False);pd.DataFrame(window_rows).to_csv(R+'/q2_lag_windows.csv',index=False);pd.DataFrame(pred_rows).to_csv(R+'/q2_lag_loho_prediction.csv',index=False)
for h,df in curves.items():df.to_csv(R+f'/curves_{h}.csv',index=False)
print('\nLAG BY HEAT\n',pd.DataFrame(lag_rows).round(3).to_string(index=False))
print('\nINFORMATIVE WINDOWS\n',pd.DataFrame(window_rows).groupby(['gas','informative']).size().to_string())
print('\nLOHO PREDICTION\n',pd.DataFrame(pred_rows).groupby(['gas','model']).test_MAE.mean().round(4).to_string())
