import argparse, os, json, time
from pathlib import Path
from functools import lru_cache
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

P=argparse.ArgumentParser();P.add_argument('--cache',type=Path,required=True);P.add_argument('--out',type=Path,required=True);P.add_argument('--heat',default='all');P.add_argument('--gas',choices=['CO2','CO','both'],default='both');P.add_argument('--intervals',default='5,15,30,60,120,180');P.add_argument('--gains',default='0.25,0.5,0.75,1.0');P.add_argument('--mode',choices=['reproduce','feedback_sensitivity'],default='reproduce');P.add_argument('--score-start-s',type=float,default=120.);P.add_argument('--prefix-s',type=float,default=60.);P.add_argument('--workers',type=int,default=4);args=P.parse_args();out=args.out;out.mkdir(exist_ok=True,parents=True)
assert args.score_start_s >= 60 + args.prefix_s
D={};hs=tuple(p.stem for p in sorted(args.cache.glob('6*.npz')));letters=dict(zip(hs,'ABCDEFGHIJ'));assert len(hs)==10

def trailing(x,n):
 x=np.asarray(x);s=np.concatenate([np.zeros_like(x[:1]),np.cumsum(x,axis=0)]);j=np.arange(1,len(x)+1);k=np.maximum(0,j-n);return (s[j]-s[k])/(j-k).reshape((-1,)+(1,)*(x.ndim-1))
def lag(x,n):return x[np.maximum(0,np.arange(len(x))-n)]
for h in hs:
 a=np.load(args.cache/(h+'.npz'));O=a['O'].astype(float);O-=O[0];t=np.arange(len(O))*.5;z=a['bands'].astype(float);c=a['contrast'].astype(float);z8=z.reshape(-1,8,4).mean(2);c8=c.reshape(-1,8,4).mean(2);br=z8.mean(1);shape=z8-br[:,None];m=t<=args.prefix_s
 x=trailing(z8,30);sh=x-x.mean(1)[:,None];q=(O-lag(O,20))/np.maximum(t-lag(t,20),.5);q60=(O-lag(O,120))/np.maximum(t-lag(t,120),.5)
 S=np.c_[x,sh,sh-sh[m].mean(0),trailing(c8,30),x-lag(x,20),x-trailing(z8,120),trailing(z8,120),trailing(z8,6)-x]
 S32=np.c_[trailing(z,30),trailing(c,30),trailing(z,30)-lag(trailing(z,30),20),S]
 pref=np.r_[x[m].mean(0).mean(),x[m].mean(1).std(),q[m].mean(),q[m].std(),sh[m].mean(0).reshape(4,2).mean(1),c[m].mean()]
 D[h]=dict(t=t,O=O,q=q,q60=q60,S=S,S32=S32,pref=pref,y=a['gas'].astype(float),br=x.mean(1),idx=np.arange(int(round(args.score_start_s*2)),len(t)))

@lru_cache(None)
def q2_stats(h,d,g):
 a=D[h];i=a['idx'][::4];t=a['t'][i];u=a['O'][i]/3000;q=a['q'][i-2*d]/5
 X=np.c_[np.ones(len(i)),t/600,u,np.maximum(u[:,None]-np.array([.15,.4,.7,.9]),0),q,q*q];y=a['y'][i,g]
 return X,y,X.T@X/len(i),X.T@y/len(i)
@lru_cache(None)
def lead(train,g):
 ds=range(6,61,6);loss=[]
 for d in ds:
  stats=[q2_stats(h,d,g) for h in train];A=sum(z[2] for z in stats);B=sum(z[3] for z in stats);pen=np.eye(A.shape[0])*.001*(len(train)-1);pen[0,0]=0
  loss.append(np.mean([np.abs(y-X@np.linalg.solve(A-aa+pen,B-bb)).mean() for X,y,aa,bb in stats]))
 return list(ds)[int(np.argmin(loss))]
@lru_cache(None)
def features(h,d):
 a=D[h];i=a['idx']-2*d;t=a['t'][i];O=a['O'][i];q=a['q'][i];u=np.maximum.accumulate(O+1.25*d*O/np.maximum(t,1));p=u/1000
 X=np.c_[p,p*p,p*p*p,np.maximum(p[:,None]-np.array([.5,1,1.5,2,2.5,3,3.5]),0)**2,q,a['q60'][i],O/np.maximum(t,1),t/100,q*p,a['q60'][i]*p]
 return u,X,a['S'][i],a['S32'][i]
@lru_cache(None)
def template(h,d,g):
 u=features(h,d)[0];y=uniform_filter1d(D[h]['y'][D[h]['idx'],g],24,mode='nearest');keep=np.r_[True,np.diff(u)>1e-7];return u[keep],y[keep]
G=np.linspace(0,5000,201)
@lru_cache(4096)
def prior(train,h,d,g):
 u=features(h,d)[0];V=np.array([np.interp(G,*template(j,d,g)) for j in train]);mu=V.mean(0);med=np.median(V,0);z=np.stack([D[j]['pref'] for j in train]);sd=np.maximum(z.std(0),1e-5);z=(z-z.mean(0))/sd;pt=np.clip((D[h]['pref']-np.stack([D[j]['pref'] for j in train]).mean(0))/sd,-3,3);s=max(mu.std(),1);coef=[]
 for v in V:
  A=np.c_[np.ones(len(mu)),mu/s];aa,bb=np.linalg.solve(A.T@A/len(mu)+np.eye(2)*.05,A.T@((v-mu)/s)/len(mu));coef.append([aa,np.log(np.clip(1+bb,.25,4))])
 aa,bb=Ridge(alpha=10).fit(z,coef).predict(pt[None])[0];sim=aa*s+np.exp(np.clip(bb,-.7,.7))*mu
 _,_,vt=np.linalg.svd(V-mu,full_matrices=False);phi=vt[:2];fp=mu+Ridge(alpha=10).fit(z,(V-mu)@phi.T).predict(pt[None])[0]@phi
 return {'median':np.interp(u,G,med),'sim':np.interp(u,G,sim),'fpcr':np.interp(u,G,fp)}
def design(h,d,mode):
 u,p,s,s32=features(h,d);phase=np.clip(u/3000,0,1.5)[:,None]
 if mode=='process':return p
 if mode=='spectral':return s
 if mode=='spectral32':return s32
 if mode=='fusion':return np.c_[p,s,s[:,:16]*(phase-.5)]
 if mode=='fusion32':return np.c_[p,s32,s[:,:16]*(phase-.5)]
 raise ValueError(mode)
def trainpredict(train,h,d,g,mode,base=None):
 X=[design(j,d,mode)[::2] for j in train];xt=design(h,d,mode);Y=[]
 for j in train:
  y=D[j]['y'][D[j]['idx'],g];Y.append((y if base is None else y-prior(tuple(k for k in train if k!=j),j,d,g)[base])[::2])
 mu=np.mean([x.mean(0) for x in X],0);sd=np.maximum(np.sqrt(np.mean([((x-mu)**2).mean(0) for x in X],0)),1e-6);X=[np.clip((x-mu)/sd,-5,5) for x in X];xt=np.clip((xt-mu)/sd,-5,5);X=[np.c_[np.ones(len(x)),x] for x in X];xt=np.c_[np.ones(len(xt)),xt]
 AA=np.array([x.T@x/len(x) for x in X]);bb=np.array([x.T@y/len(y) for x,y in zip(X,Y)]);reg=np.eye(xt.shape[1]);reg[0,0]=1e-8;outp={};a0=0 if base is None else prior(train,h,d,g)[base]
 for alpha in (.01,.03,.1,.3,1.,3.,10.):
  coef=np.linalg.solve(AA.mean(0)+alpha*reg,bb.mean(0));v=xt@coef
  for strength in ((1.,) if base is None else (.25,.5,.75,1.,1.25)):
   outp[f'{mode}_{base or "direct"}_ridge{alpha}_s{strength}']=a0+strength*v
  if base is None and mode in ('spectral','fusion'):
   co=np.linalg.solve(AA+alpha*reg,bb[...,None]).squeeze(-1);outp[f'{mode}_experts_ridge{alpha}']=xt@co.mean(0)
 if base is None and mode in ('spectral','fusion'):
  xx=np.concatenate([x[:,1:] for x in X]);yy=np.concatenate(Y)
  for n in (3,8):
   mdl=PLSRegression(n_components=n,scale=False,max_iter=100).fit(xx,yy);outp[f'{mode}_PLS{n}']=mdl.predict(xt[:,1:]).ravel()
 return outp
@lru_cache(128)
def candidates(train,h,g):
 assert h not in train
 d=lead(train,g);pr=prior(train,h,d,g);p={'process_template':pr['median'],'prior_sim':pr['sim'],'prior_fpcr':pr['fpcr']}
 for mode in ('process','spectral','spectral32','fusion','fusion32'):p.update(trainpredict(train,h,d,g,mode))
 for mode in ('process','fusion','fusion32'):
  for base in (('median',) if mode=='process' else ('median','sim' if g==0 else 'fpcr')):p.update(trainpredict(train,h,d,g,mode,base))
 return {k:np.maximum(v,0) for k,v in p.items()},d
def family(n):
 if n.startswith('prior'):return 'prior'
 if n.startswith('process'):return 'process'
 if n.startswith('spectral'):return 'spectral'
 if n.startswith('fusion'):return 'fusion'
 raise ValueError(n)
def score(y,v):
 e=y-v;return dict(r2=1-np.sum(e**2)/np.sum((y-y.mean())**2),rmse=np.sqrt(np.mean(e**2)),mae=np.mean(abs(e)))
from concurrent.futures import ProcessPoolExecutor
COARSE=np.arange(0,1.5001,.25)
FINE=np.arange(0,1.5001,.05)
STEPS=np.arange(1,31)

def indices(n,d,step):
 k=np.arange(n); pub=k-2*d; last=(pub//step)*step; ok=last>=0
 assert np.all(last[ok]<=pub[ok]);return last,ok

def feedback(v,y,d,step,gain):
 last,ok=indices(len(y),d,step);r=np.zeros(len(y));r[ok]=y[last[ok]]-v[last[ok]]
 return np.maximum(0,v+gain*r)

def arrays(train,h,g):
 p,d=candidates(train,h,g);y=D[h]['y'][D[h]['idx'],g]
 names=sorted(n for n in p if family(n) in ('fusion','process','prior'))
 return names,np.stack([p[n] for n in names]),y,d

def grid(v,y,d,steps,gains):
 den=np.sum((y-y.mean())**2); out=np.empty((len(steps),len(v),len(gains)))
 for a,step in enumerate(steps):
  last,ok=indices(len(y),d,step);r=np.zeros_like(v);r[:,ok]=y[last[ok]][None,:]-v[:,last[ok]]
  e=np.maximum(0,v[:,None,:]+gains[None,:,None]*r[:,None,:])-y
  out[a]=1-np.sum(e*e,axis=2)/den
 return out

def fitfold(task):
 h,g=task; tr=tuple(j for j in hs if j!=h); data=[]
 with threadpool_limits(limits=1):
  for j in tr:
   nn,v,y,d=arrays(tuple(k for k in tr if k!=j),j,g);
   if data: assert nn==names
   names=nn;data.append((v,y,d))
  coarse=np.mean([grid(v,y,d,STEPS,COARSE) for v,y,d in data],0)
  fam=np.array([family(n) for n in names]); ix_f=np.flatnonzero(fam=='fusion');ix_p=np.flatnonzero(fam=='process');ix_prior=np.flatnonzero(fam=='prior')
  testnames,test,y,d=arrays(tr,h,g);assert testnames==names
  rows=[]
  for label,ids in [('fusion',ix_f),('process',ix_p),('prefix_prior',ix_prior)]:
   b=ids[np.argmax(coarse[0,ids,0])];v=test[b];rows.append(dict(heat=h,gas=('CO2','CO')[g],family=label,interval_s=0.,gain=0.,model=names[b],lead_s=d,inner_r2=coarse[0,b,0],n=len(y),**score(y,v)))
   for k,step in enumerate(STEPS):
    chosen=ids[np.argsort(coarse[k,ids].max(axis=1))[-3:]]
    fine=np.mean([grid(v[chosen],yy,dd,[step],FINE)[0] for v,yy,dd in data],0)
    ii,jj=np.unravel_index(np.argmax(fine),fine.shape);b=chosen[ii];gain=FINE[jj]
    vv=feedback(test[b],y,d,step,gain)
    rows.append(dict(heat=h,gas=('CO2','CO')[g],family=label,interval_s=step/2,gain=gain,model=names[b],lead_s=d,inner_r2=fine[ii,jj],n=len(y),**score(y,vv)))
  fallback=np.mean([D[j]['y'][D[j]['idx'],g].mean() for j in tr])
  for step in STEPS:
   last,ok=indices(len(y),d,step);v=np.full(len(y),fallback);v[ok]=y[last[ok]]
   rows.append(dict(heat=h,gas=('CO2','CO')[g],family='persistence',interval_s=step/2,gain=1.,model='last_observation',lead_s=d,inner_r2=np.nan,n=len(y),**score(y,v)))
  for label in ['fusion','process','prefix_prior']:
   opts=[r for r in rows if r['family']==label and r['interval_s']>0];best=max(opts,key=lambda r:r['inner_r2']);rows.append(dict(best,family=label+'_nested_interval'))
 return rows

def main():
 tasks=[(h,g) for h in hs for g in (0,1)]; rows=[]
 with ProcessPoolExecutor(max_workers=args.workers) as ex:
  for task,rr in zip(tasks,ex.map(fitfold,tasks)):
   rows+=rr;pd.DataFrame(rows).to_csv(out/'fold_metrics_partial.csv',index=False);print('DONE',task,flush=True)
 df=pd.DataFrame(rows);df['score_start_s']=args.score_start_s;df['prefix_s']=args.prefix_s;df.to_csv(out/'fold_metrics.csv',index=False)
 regular=df[~df.family.str.endswith('_nested_interval')];s=regular.groupby(['gas','family','interval_s']).agg(heats=('heat','nunique'),macro_r2=('r2','mean'),macro_rmse=('rmse','mean'),worst_r2=('r2','min'),median_gain=('gain','median')).reset_index();s.to_csv(out/'summary.csv',index=False)
 nested=df[df.family.str.endswith('_nested_interval')];print(nested.groupby(['gas','family']).agg(r2=('r2','mean'),rmse=('rmse','mean'),worst=('r2','min'),interval=('interval_s','median')).to_string(),flush=True)
 print(s[s.family=='fusion'].sort_values(['gas','macro_r2'],ascending=[True,False]).groupby('gas').head(5).to_string(index=False),flush=True)
if __name__=='__main__':main()
