from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter1d
from sklearn.preprocessing import SplineTransformer
from threadpoolctl import threadpool_limits
from functools import lru_cache
import argparse
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path(__file__).parent.parent/'q2_emission_rerun/test_process.csv');p.add_argument('--out',type=Path,default=Path(__file__).parent);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
f=pd.read_csv(args.input);D={h:(a.O.to_numpy(),a[['CO2','CO']].to_numpy()) for h,a in f.sort_values(['attachment','row']).groupby('heat',sort=False)};H=tuple(D);assert len(H)==10 and sum(len(a) for a,b in D.values())==13589 and H[4]=='6C02057'
B=[([0,0,-500,1,-500,1],[120,120,4000,3000,4500,3000]),([0,-500,1,0,1],[150,4000,3000,5500,3000])]; P=[[27,14,124,160,950,400],[45,833,180,3294,120]]
def fun(g,x,p):
 if g==0:
  a,b,v,w,u,s=p;return a*expit((x-v)/w)-b*expit((x-u)/s)
 c,v,w,u,s=p;return c*expit((x-v)/w)*expit((u-x)/s)
def score(y,z):
 e=y-z;return 1-(e@e)/np.sum((y-y.mean())**2),np.sqrt(np.mean(e**2)),np.mean(abs(e))
@lru_cache(None)
def sigmoid(train,g):
 x=np.concatenate([D[h][0][120::6] for h in train]);y=np.concatenate([D[h][1][120::6,g] for h in train]);w=np.concatenate([np.full(len(D[h][0][120::6]),1/np.sqrt(len(D[h][0][120::6]))) for h in train]);
 def res(p):
  z=fun(g,x,p);return np.r_[(z-y)*w,np.minimum(z,0)*w*3]
 return least_squares(res,P[g],bounds=B[g],max_nfev=100,ftol=1e-6).x
@lru_cache(None)
def feat(h,sc):
 o,y=D[h];dt=1 if sc=='H1' and h==H[4] else .5;t=np.arange(len(o))*dt;u=(o-o[0])/1000
 qs=[]
 for sec in (10,30,60):
  j=np.maximum(0,np.arange(len(o))-int(sec/dt));qs.append((o-o[j])/np.maximum(t-t[j],dt)/5)
 a=np.c_[t/600,u,np.clip(u/np.maximum(t,dt)*200,0,3),*qs]
 return u,a,u/max(u[-1],.01)
@lru_cache(None)
def curve(h,g,kind,band):
 o,y=D[h];x=feat(h,'H1')[0 if kind=='absolute' else 2][120:];z=y[120:,g];x=np.maximum.accumulate(x);k=np.r_[np.diff(x)>0,True];x=x[k];z=z[k];grid=np.arange(0,4.501,.005) if kind=='absolute' else np.linspace(0,1,401);v=np.interp(grid,x,z);v=gaussian_filter1d(v,band/(5 if kind=='absolute' else 8.5),mode='nearest') if band else v;return grid,v,(grid>=x[0])&(grid<=x[-1])
@lru_cache(None)
def design(h,sc,kind,n):
 u,a,r=feat(h,sc);x=(r if kind=='relative' else u)[:,None];kn=np.linspace(0,1 if kind=='relative' else 4,n)[:,None];b=SplineTransformer(knots=kn,degree=3,include_bias=True,extrapolation='constant').fit_transform(x)
 if kind in ('absolute','relative'):return b
 z=a[:,[0,2,3,4,5]];return np.c_[b,z,(b[:,:,None]*z[:,None,:]).reshape(len(b),-1)]
@lru_cache(None)
def gram(h,g,sc,kind,n):
 x=design(h,sc,kind,n)[120::3];y=D[h][1][120::3,g];return x.T@x/len(x),x.T@y/len(x)
def pred_all(train,tests,g,sc):
 out={h:{} for h in tests};p=sigmoid(train,g)
 for h in tests:out[h]['sigmoid']=fun(g,D[h][0],p)
 for kind in ('absolute','relative'):
  for band in (0,40,100,200):
   cv=[curve(h,g,kind,band) for h in train];grid=cv[0][0];V=np.array([v for _,v,_ in cv]);mask=np.array([m for _,_,m in cv]);vm=np.where(mask,V,np.nan)
   supported=mask.sum(0)>0;mean=V.mean(0);med=np.median(V,axis=0);mean[supported]=np.nanmean(vm[:,supported],axis=0);med[supported]=np.nanmedian(vm[:,supported],axis=0)
   for h in tests:
    x=feat(h,sc)[0 if kind=='absolute' else 2];out[h][f'{kind}_mean_{band}']=np.interp(x,grid,mean);out[h][f'{kind}_median_{band}']=np.interp(x,grid,med)
 for kind in ('absolute','relative','process'):
  for n in (5,9,15):
   a,b=zip(*(gram(h,g,sc,kind,n) for h in train));A=np.mean(a,axis=0);B=np.mean(b,axis=0)
   for alpha in (.0003,.003,.03,.3):
    c=np.linalg.solve(A+alpha*np.eye(len(B)),B)
    for h in tests:out[h][f'{kind}_spline_{n}_{alpha}']=design(h,sc,kind,n)@c
 return {h:{k:np.maximum(v,0) for k,v in z.items()} for h,z in out.items()}
def family(k):return 'absolute' if k=='sigmoid' else k.split('_')[0]
metrics=[];predictions=[];inners=[]
with threadpool_limits(limits=1):
 for sc in ('H1','H0'):
  for g,gas in enumerate(('CO2','CO')):
   cache={};outer={}
   for i,h in enumerate(H):
    outer[h]=pred_all(tuple(k for k in H if k!=h),(h,),g,sc)[h]
    for j in range(i):
     v=H[j];p=pred_all(tuple(k for k in H if k not in (h,v)),(h,v),g,sc);cache[h,v]=p[v];cache[v,h]=p[h]
   for h in H:
    y=D[h][1][120:,g];p=outer[h];rank={k:np.mean([score(D[v][1][120:,g],cache[h,v][k][120:])[0] for v in H if v!=h]) for k in p}
    chosen={'sigmoid':'sigmoid'}
    for fam in ('absolute','relative','process'):chosen[fam]=max((k for k in p if family(k)==fam),key=rank.get)
    for k,v in rank.items():inners.append((sc,gas,h,k,v))
    for fam,k in chosen.items():
     z=p[k][120:];r,rm,ma=score(y,z);metrics.append((sc,gas,h,fam,k,rank[k],r,rm,ma));predictions.extend(zip([sc]*len(y),[gas]*len(y),[h]*len(y),[fam]*len(y),D[h][0][120:],y,z))
   print(sc,gas,pd.DataFrame(metrics,columns=['scenario','gas','heat','family','model','inner','r2','rmse','mae']).query('scenario==@sc and gas==@gas').groupby('family').r2.mean().to_dict(),flush=True)
cols=['scenario','gas','heat','family','model','inner','r2','rmse','mae'];m=pd.DataFrame(metrics,columns=cols);m.to_csv(args.out/'metrics.csv',index=False);pd.DataFrame(inners,columns=['scenario','gas','heat','model','inner_r2']).to_csv(args.out/'inner.csv',index=False);p=pd.DataFrame(predictions,columns=['scenario','gas','heat','family','oxygen','observed','predicted']);p.to_csv(args.out/'predictions.csv.gz',index=False);s=m.groupby(['scenario','gas','family'])[['r2','rmse','mae']].mean().reset_index();s['pooled_r2']=[score(v.observed.to_numpy(),v.predicted.to_numpy())[0] for _,v in p.groupby(['scenario','gas','family'])];s.to_csv(args.out/'summary.csv',index=False);print(s.to_string(index=False))
