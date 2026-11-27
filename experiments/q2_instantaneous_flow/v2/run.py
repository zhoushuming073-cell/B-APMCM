from pathlib import Path
from functools import lru_cache
import numpy as np,pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
import argparse
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path(__file__).parent/'process.csv');p.add_argument('--out',type=Path,default=Path(__file__).parent);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
f=pd.read_csv(a.input);D={h:(d.O.to_numpy(),d[['CO2','CO']].to_numpy()) for h,d in f.sort_values(['attachment','row']).groupby('heat',sort=False)};H=tuple(D);assert len(H)==10 and sum(len(x) for x,y in D.values())==13589 and H[4]=='6C02057'
G=np.arange(0,4001,10.);Y={};X={};Z={};M={}
for h,(o,y) in D.items():
 x=np.maximum.accumulate(o-o[0]);X[h]=x;idx=np.r_[np.diff(x)>0,True];Y[h]=np.column_stack([gaussian_filter1d(np.interp(G,x[idx],y[idx,g]),6,mode='nearest') for g in range(2)]);M[h]=(G>=x[120])&(G<=x[-1]);
 for sc in ('H0','H1'):
  dt=1 if sc=='H1' and h==H[4] else .5;t=np.arange(len(o))*dt;j=120;q=(x[20:j+1]-x[:j-19])/(20*dt);q=np.maximum(q,0)
  inv=np.interp(np.linspace(0,3000,16),x[idx],t[idx]);rates=200/np.maximum(np.diff(inv),1)
  Z[h,sc,'prefix']=np.r_[x[j]/t[j],np.mean(q),np.std(q),np.quantile(q,[.1,.5,.9]),np.diff(x[np.arange(0,121,20)])/(20*dt)]
  Z[h,sc,'full']=np.r_[Z[h,sc,'prefix'],x[-1],t[-1],rates,np.quantile(rates,[.1,.5,.9])]
def score(h,g,p):
 y=D[h][1][120:,g];e=y-np.maximum(p,0);return 1-e@e/np.sum((y-y.mean())**2),np.sqrt(np.mean(e**2)),np.mean(abs(e))
@lru_cache(None)
def mu(train):
 V=np.array([Y[h] for h in train]);m=np.array([M[h] for h in train]);n=m.sum(0);z=V.mean(0);ok=n>0;z[ok]=(V*m[:,:,None]).sum(0)[ok]/n[ok,None];return z
@lru_cache(None)
def stats(train,sc,kind):
 z=np.stack([Z[h,sc,kind] for h in train]);c=z.mean(0);s=np.maximum(z.std(0),1e-5);return (z-c)/s,c,s
@lru_cache(None)
def align(train,g):
 b=mu(train)[:,g];scale=max(np.std(b),1);params=[]
 for h in train:
  x=X[h][120::4];y=D[h][1][120::4,g]
  def residual(p):
   aa,ss,c,d=p;pr=np.exp(aa)*np.interp(np.exp(ss)*x+c*500,G,b)+d*scale;return np.r_[(pr-y)/scale/np.sqrt(len(y)),.12*p]
  v=least_squares(residual,[0.,0.,0.,0.],bounds=([-1,-.35,-1.2,-1],[1,.35,1.2,1]),max_nfev=60).x;params.append(v)
 return b,np.array(params)
@lru_cache(None)
def pred(train,h,g,sc):
 assert h not in train;x=X[h][120:];V=np.array([Y[j][:,g] for j in train]);b=mu(train)[:,g];out={'shared':np.interp(x,G,b)}
 for kind in ('prefix','full'):
  z,c,s=stats(train,sc,kind);zt=np.clip((Z[h,sc,kind]-c)/s,-3,3);dist=np.mean((z-zt)**2,1)
  for temp in (.25,1,4):
   w=np.exp(-(dist-dist.min())/temp);w/=w.sum();m=np.array([M[j] for j in train]);ww=w[:,None]*m;v=(ww*V).sum(0)/np.maximum(ww.sum(0),1e-12);v=np.where(ww.sum(0)>1e-10,v,b)
   for shrink in (.25,.75):out[f'{kind}_kernel_{temp}_{shrink}']=np.interp(x,G,b+shrink*(v-b))
  for K in (2,3):
   lab=KMeans(n_clusters=K,n_init=5,random_state=5).fit_predict(V);dd=np.array([np.mean((z[lab==k].mean(0)-zt)**2) for k in range(K)]);w=np.exp(-(dd-dd.min()));w/=w.sum();v=sum(w[k]*V[lab==k].mean(0) for k in range(K));out[f'{kind}_group_{K}']=np.interp(x,G,.5*b+.5*v)
  for joint in (False,True):
   A=np.stack([Y[j].T.ravel() for j in train]) if joint else V;avg=A.mean(0);_,sig,vt=np.linalg.svd(A-avg,full_matrices=False)
   for k in (1,2,3):
    for alpha in (3,15,60):
     coeff=Ridge(alpha=alpha).fit(z,(A-avg)@vt[:k].T).predict(zt[None])[0];v=np.atleast_1d(coeff)@vt[:k];v=v.reshape(2,-1)[g] if joint else v;out[f'{kind}_factor_{int(joint)}_{k}_{alpha}']=np.interp(x,G,b+v)
  base,pars=align(train,g)
  for alpha in (3,15,60):
   v=Ridge(alpha=alpha).fit(z,pars).predict(zt[None])[0];v=np.clip(v,[-1,-.35,-1.2,-1],[1,.35,1.2,1]);aa,ss,cc,dd=v;out[f'{kind}_SIM_{alpha}']=np.exp(aa)*np.interp(np.exp(ss)*x+cc*500,G,base)+dd*max(base.std(),1)
 return {k:np.maximum(v,0) for k,v in out.items()}
def fam(n):return 'shared' if n=='shared' else '_'.join(n.split('_')[:2])
metrics=[];inner=[];preds=[];oracles=[]
with threadpool_limits(limits=1):
 for sc in ('H1','H0'):
  for g,gas in enumerate(('CO2','CO')):
   for h in H:
    train=tuple(j for j in H if j!=h);out=pred(train,h,g,sc);r={k:[] for k in out}
    for v in train:
     pp=pred(tuple(j for j in train if j!=v),v,g,sc)
     for k,pr in pp.items():r[k].append(score(v,g,pr)[0])
    ranks={k:np.mean(v) for k,v in r.items()};inner.extend((sc,gas,h,k,v) for k,v in ranks.items());choices={family:max((k for k in out if fam(k)==family),key=ranks.get) for family in sorted(set(map(fam,out)))}
    for kind in ('prefix','full'):choices[kind+'_selected']=max((k for k in out if k.startswith(kind) or k=='shared'),key=ranks.get)
    y=D[h][1][120:,g];x=X[h][120:];b=mu(train)[:,g];scale=max(b.std(),1)
    def res(p):
     aa,ss,cc,dd=p;return np.exp(aa)*np.interp(np.exp(ss)*x+cc*500,G,b)+dd*scale-y
    v=least_squares(res,[0,0,0,0],bounds=([-1,-.35,-1.2,-1],[1,.35,1.2,1]),max_nfev=100).x;oracles.append((sc,gas,h,*v,score(h,g,y+res(v))[0]))
    for family,k in sorted(choices.items()):
     pp=out[k];rr,rm,ma=score(h,g,pp);metrics.append((sc,gas,h,family,k,ranks[k],rr,rm,ma));preds.extend(zip([sc]*len(x),[gas]*len(x),[h]*len(x),[family]*len(x),x,y,pp))
   print(sc,gas,flush=True)
m=pd.DataFrame(metrics,columns=['scenario','gas','heat','family','model','inner_r2','r2','rmse','mae']);m.to_csv(a.out/'metrics.csv',index=False);pd.DataFrame(inner,columns=['scenario','gas','heat','model','inner_r2']).to_csv(a.out/'inner.csv',index=False);pr=pd.DataFrame(preds,columns=['scenario','gas','heat','family','oxygen_net','observed','predicted']);pr.to_csv(a.out/'predictions.csv.gz',index=False);s=m.groupby(['scenario','gas','family'])[['r2','rmse','mae']].mean().reset_index();s.to_csv(a.out/'summary.csv',index=False);pd.DataFrame(oracles,columns=['scenario','gas','heat','log_amplitude','log_stretch','shift_500','offset_scaled','r2']).to_csv(a.out/'oracle.csv',index=False);print(s.to_string(index=False))
