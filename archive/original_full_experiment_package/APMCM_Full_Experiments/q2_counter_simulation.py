import os,glob
import numpy as np,pandas as pd
from scipy.interpolate import BSpline,PchipInterpolator
from scipy.optimize import lsq_linear
R='/mnt/data/APMCM_Full_Experiments';rows=[]
def fitflow(O,t,spacing=10,lam=10):
 knots=np.r_[[t[0]]*4,np.arange(t[0]+spacing,t[-1],spacing),[t[-1]]*4];p=len(knots)-4;bs=BSpline(knots,np.eye(p),3);B=bs(t);anti=bs.antiderivative();A=anti(t)-anti(t[0]);D2=np.diff(np.eye(p),n=2,axis=0);w=np.ones(len(t))
 for _ in range(2):
  M=np.vstack([A*np.sqrt(w)[:,None],np.sqrt(lam)*D2]);y=np.r_[O-O[0],np.zeros(len(D2))];f=lsq_linear(M,y,bounds=(0,np.inf),max_iter=80);c=f.x;r=O-O[0]-A@c;s=max(1.4826*np.median(abs(r-np.median(r))),.1);w=np.minimum(1,1.345*s/np.maximum(abs(r),1e-9))
 return B@c

def slope(O,t,win,causal=False):
 n=len(O);k=int(round(win/.5));q=np.full(n,np.nan)
 for i in range(n):
  a=max(0,i-k+1) if causal else max(0,i-k//2);b=i+1 if causal else min(n,i+k//2+1)
  if b-a<3:continue
  x=t[a:b];y=O[a:b];xc=x-x.mean();q[i]=max(0,xc@(y-y.mean())/(xc@xc))
 return q
for f in sorted(glob.glob('/mnt/data/heats_*.npz')):
 z=np.load(f);h=os.path.basename(f)[6:-4];n=len(z['seq']);t=np.arange(n)*.5;O=z['proc'][:,0].astype(float);updates=np.r_[0,np.flatnonzero(np.diff(O)>0)+1]
 rg=np.random.default_rng(20261201+sum(map(ord,h)))
 scenarios={'step':np.where(t<.15*t[-1],3.0,np.where(t<.45*t[-1],5.2,np.where(t<.72*t[-1],4.1,4.8))),
 'pulse':4.7+2.0*np.exp(-((t-.42*t[-1])/2.5)**2)-1.5*np.exp(-((t-.68*t[-1])/3.0)**2),
 'oscillatory':4.7+.7*np.sin(2*np.pi*t/45)+.25*np.sin(2*np.pi*t/8)}
 for scen,qtrue in scenarios.items():
  otr=np.r_[0,np.cumsum((qtrue[:-1]+qtrue[1:])*.25)]
  # sample/hold at the empirical update times and retain the empirical update interval structure
  obs=np.empty(n);vals=np.round(otr[updates],2)
  for k,ix in enumerate(updates):
   end=updates[k+1] if k+1<len(updates) else n;obs[ix:end]=vals[k]
  # methods
  qint=np.zeros(n)
  for k in range(len(updates)-1):
   a,b=updates[k],updates[k+1];v=(obs[b]-obs[a])/(t[b]-t[a]);qint[a:b]=max(0,v)
  qint[updates[-1]:]=qint[max(0,updates[-1]-1)] if updates[-1]>0 else 0
  pred={'update_interval_slope':qint,'PCHIP_derivative':np.maximum(0,PchipInterpolator(t,obs).derivative()(t)),'nonnegative_integral_spline':fitflow(obs,t,10,10),'centered_slope_2s':slope(obs,t,2),'centered_slope_5s':slope(obs,t,5),'centered_slope_10s':slope(obs,t,10),'trailing_slope_10s':slope(obs,t,10,True)}
  valid=(t>=10)&(t<=t[-1]-10);rapid=valid&(np.abs(np.gradient(qtrue,.5))>.08)
  for method,qhat in pred.items():
   v=valid&np.isfinite(qhat);r=rapid&np.isfinite(qhat)
   rows.append({'heat':h,'scenario':scen,'method':method,'MAE':float(np.mean(abs(qhat[v]-qtrue[v]))),'RMSE':float(np.sqrt(np.mean((qhat[v]-qtrue[v])**2))),'rapid_MAE':float(np.mean(abs(qhat[r]-qtrue[r]))) if r.sum() else np.nan,'peak_error':float(abs(np.max(qhat[v])-np.max(qtrue[v]))),'n':int(v.sum())})
 print('simulated update mechanism',h,'updates',len(updates),flush=True)
df=pd.DataFrame(rows);df.to_csv(R+'/q2_counter_simulation_allheats.csv',index=False)
print(df.groupby(['scenario','method'])[['MAE','rapid_MAE','peak_error']].mean().round(3).to_string())
