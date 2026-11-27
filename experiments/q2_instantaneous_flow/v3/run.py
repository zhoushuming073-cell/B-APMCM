from pathlib import Path
from functools import lru_cache
import argparse
import numpy as np,pandas as pd
from scipy.ndimage import gaussian_filter1d as smooth
from threadpoolctl import threadpool_limits
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path(__file__).parent/'process.csv');p.add_argument('--out',type=Path,default=Path(__file__).parent);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
f=pd.read_csv(a.input);D={h:(d.O.to_numpy(),d[['CO2','CO']].to_numpy()) for h,d in f.sort_values(['attachment','row']).groupby('heat',sort=False)};H=tuple(D);assert len(H)==10 and sum(len(x) for x,y in D.values())==13589 and H[4]=='6C02057'
G=np.arange(0,4001,10.);Y={};X={};M={};Z={}
def interp(v,x):return np.column_stack([np.interp(x,G,v[:,i]) for i in range(v.shape[1])])
for h,(o,y) in D.items():
 x=np.maximum.accumulate(o-o[0]);X[h]=x;k=np.r_[np.diff(x)>0,True];Y[h]=np.column_stack([smooth(np.interp(G,x[k],y[k,g]),6,mode='nearest') for g in range(2)]);M[h]=(G>=x[120])&(G<=x[-1]);k=np.r_[True,np.diff(x)>0]
 for sc in ('H0','H1'):
  dt=1 if sc=='H1' and h==H[4] else .5;t=np.arange(len(o))*dt;j=np.arange(len(o));q=[]
  for sec in (10,30,60,120):
   l=np.maximum(0,j-int(sec/dt));q.append((x-x[l])/np.maximum(t-t[l],dt))
  q=np.array(q).T;avg=x/np.maximum(t,dt);pref=np.r_[avg[120],q[120,1:],np.std(q[:121,0])];local=np.c_[q,avg,q[:,0]-q[:,2],q[:,1]-q[:,3],t/120];pref=np.tile(pref,(len(x),1))
  for kind,v in [('local',local),('prefix',pref),('combined',np.c_[local,pref])]:Z[h,sc,kind]=(np.column_stack([np.interp(G,x[k],v[k,i]) for i in range(v.shape[1])]),v[120:])
def score(h,g,p):
 y=D[h][1][120:,g];e=y-p;return 1-e@e/np.sum((y-y.mean())**2),np.sqrt(np.mean(e**2)),np.mean(abs(e))
@lru_cache(None)
def base(tr,g):
 V=np.array([Y[h][:,g] for h in tr]);m=np.array([M[h] for h in tr]);n=m.sum(0);b=V.mean(0);ok=n>0;b[ok]=(V*m).sum(0)[ok]/n[ok];return V,m,b
@lru_cache(None)
def fit(tr,g,sc,kind):
 V,m,b=base(tr,g);z=np.stack([Z[h,sc,kind][0] for h in tr]);n=m.sum(0);w=m/np.maximum(n,1);mu=np.einsum('hn,hnp->np',w,z);sd=np.sqrt(np.maximum(np.einsum('hn,hnp->np',w,(z-mu)**2),.01));zs=np.where(m[:,:,None],np.clip((z-mu)/sd,-4,4),0);res=(V-b)*m;_,_,vt=np.linalg.svd(res,full_matrices=False);gram=np.einsum('hnp,hnq->npq',zs,zs)/np.maximum(n[:,None,None],1);rhs=np.einsum('hnp,hn->np',zs,res)/np.maximum(n[:,None],1);co={}
 for bw in (100.,300.):
  A=smooth(gram,bw/10,axis=0,mode='nearest');B=smooth(rhs,bw/10,axis=0,mode='nearest')
  for alpha in (.3,3.,15.):co[f'varying_{kind}_{bw}_{alpha}']=np.linalg.solve(A+alpha*np.eye(A.shape[-1]),B[...,None]).squeeze(-1)
 for K in (1,2):
  phi=vt[:K].T
  for bw in (200.,500.):
   PP=np.einsum('nk,nl->nkl',phi,phi);A=smooth(m[:,:,None,None]*PP[None],bw/10,axis=1,mode='nearest');B=smooth(res[:,:,None]*phi[None],bw/10,axis=1,mode='nearest');load=np.linalg.solve(A+.0005*np.eye(K),B[...,None]).squeeze(-1)
   load=(load-(load*m[:,:,None]).sum(0)[None]/np.maximum(n[None,:,None],1))*m[:,:,None]
   rhs=np.einsum('hnp,hnk->npk',zs,load)/np.maximum(n[:,None,None],1);AA=smooth(gram,bw/10,axis=0,mode='nearest');BB=smooth(rhs,bw/10,axis=0,mode='nearest')
   for alpha in (1.,10.):
    cc=np.linalg.solve(AA+alpha*np.eye(AA.shape[-1]),BB);co[f'factor_{kind}_{K}_{bw}_{alpha}']=np.einsum('npk,nk->np',cc,phi)
 return V,m,b,mu,sd,zs,n,co
@lru_cache(None)
def predict(tr,h,g,sc):
 assert h not in tr;x=X[h][120:];V,m,b=base(tr,g);out={'shared':np.interp(x,G,b)}
 for kind in ('local','prefix','combined'):
  V,m,b,mu,sd,zs,n,co=fit(tr,g,sc,kind);zt=np.clip((Z[h,sc,kind][1]-interp(mu,x))/interp(sd,x),-4,4);valid=np.interp(x,G,(n>=3).astype(float))>=.999
  ztrain=np.stack([interp(z,x) for z in zs]);vtrain=np.stack([np.interp(x,G,v) for v in V]);mask=np.stack([np.interp(x,G,mm.astype(float))>=.999 for mm in m]);dist=((ztrain-zt[None])**2).mean(2)
  for temp in (.5,2.):
   w=np.exp(-dist/temp)*mask;den=w.sum(0);v=(w*vtrain).sum(0)/np.maximum(den,1e-12);bt=out['shared'];v=np.where((den>1e-8)&valid,v,bt)
   for strength in (.25,.75):out[f'kernel_{kind}_{temp}_{strength}']=bt+strength*(v-bt)
  for k,c in co.items():out[k]=out['shared']+np.where(valid,np.sum(interp(c,x)*zt,axis=1),0)
 return {k:np.maximum(v,0) for k,v in out.items()}
def family(k):return k.split('_')[0]
metrics=[];inners=[];preds=[]
with threadpool_limits(limits=1):
 for sc in ('H1','H0'):
  for g,gas in enumerate(('CO2','CO')):
   for h in H:
    tr=tuple(j for j in H if j!=h);out=predict(tr,h,g,sc);r={k:[] for k in out}
    for v in tr:
     pp=predict(tuple(j for j in tr if j!=v),v,g,sc)
     for k,pr in pp.items():r[k].append(score(v,g,pr)[0])
    ranks={k:np.mean(v) for k,v in r.items()};inners.extend((sc,gas,h,k,v) for k,v in ranks.items());choices={fam:max((k for k in out if family(k)==fam),key=ranks.get) for fam in sorted(set(map(family,out)))};choices['selected']=max(ranks,key=ranks.get)
    for fam,k in choices.items():
     pr=out[k];rr,rm,ma=score(h,g,pr);metrics.append((sc,gas,h,fam,k,ranks[k],rr,rm,ma));preds.extend(zip([sc]*len(pr),[gas]*len(pr),[h]*len(pr),[fam]*len(pr),X[h][120:],D[h][1][120:,g],pr))
   print(sc,gas,pd.DataFrame(metrics,columns=['scenario','gas','heat','family','model','inner_r2','r2','rmse','mae']).query('scenario==@sc and gas==@gas').groupby('family').r2.mean().to_dict(),flush=True);predict.cache_clear();fit.cache_clear()
m=pd.DataFrame(metrics,columns=['scenario','gas','heat','family','model','inner_r2','r2','rmse','mae']);m.to_csv(a.out/'metrics.csv',index=False);pd.DataFrame(inners,columns=['scenario','gas','heat','model','inner_r2']).to_csv(a.out/'inner.csv',index=False);pr=pd.DataFrame(preds,columns=['scenario','gas','heat','family','oxygen_net','observed','predicted']);pr.to_csv(a.out/'predictions.csv.gz',index=False);s=m.groupby(['scenario','gas','family'])[['r2','rmse','mae']].mean().reset_index();s['pooled_r2']=[1-((v.observed-v.predicted)**2).sum()/((v.observed-v.observed.mean())**2).sum() for _,v in pr.groupby(['scenario','gas','family'])];s.to_csv(a.out/'summary.csv',index=False);print(s.to_string(index=False))
