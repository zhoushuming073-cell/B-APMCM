"""Q3 causal spectral soft sensor, nested LOHO, H0. Run with --cache prior exact NPZ cache."""
import argparse, os
from pathlib import Path
from functools import lru_cache
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy.special import softmax
from threadpoolctl import threadpool_limits
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import Ridge
P=argparse.ArgumentParser();P.add_argument('--cache',type=Path,required=True);P.add_argument('--out',type=Path,required=True);P.add_argument('--lag-scheme',choices=['spectral_lag','q2_fixed'],default='spectral_lag');args=P.parse_args();out=args.out;out.mkdir(exist_ok=True,parents=True)
D={};hs=tuple(p.stem for p in sorted(args.cache.glob('6*.npz')));letters=dict(zip(hs,'ABCDEFGHIJ'));assert len(hs)==10

def trailing(x,n):
 x=np.asarray(x);s=np.concatenate([np.zeros_like(x[:1]),np.cumsum(x,axis=0)]);j=np.arange(1,len(x)+1);k=np.maximum(0,j-n);return (s[j]-s[k])/(j-k).reshape((-1,)+(1,)*(x.ndim-1))
def lag(x,n):return x[np.maximum(0,np.arange(len(x))-n)]
for h in hs:
 a=np.load(args.cache/(h+'.npz'));O=a['O'].astype(float);O-=O[0];t=np.arange(len(O))*.5;z=a['bands'].astype(float);c=a['contrast'].astype(float);z8=z.reshape(-1,8,4).mean(2);c8=c.reshape(-1,8,4).mean(2);br=z8.mean(1);shape=z8-br[:,None];m=t<=60
 x=trailing(z8,30);sh=x-x.mean(1)[:,None];q=(O-lag(O,20))/np.maximum(t-lag(t,20),.5);q60=(O-lag(O,120))/np.maximum(t-lag(t,120),.5)
 # spectral only: no O or time covariate; all windows trailing, first 60 s available by first publish.
 S=np.c_[x,sh,sh-sh[m].mean(0),trailing(c8,30),x-lag(x,20),x-trailing(z8,120),trailing(z8,120),trailing(z8,6)-x]
 S32=np.c_[trailing(z,30),trailing(c,30),trailing(z,30)-lag(trailing(z,30),20),S]
 pref=np.r_[x[m].mean(0).mean(),x[m].mean(1).std(),q[m].mean(),q[m].std(),sh[m].mean(0).reshape(4,2).mean(1),c[m].mean()]
 D[h]=dict(t=t,O=O,q=q,q60=q60,S=S,S32=S32,pref=pref,y=a['gas'].astype(float),br=x.mean(1),idx=np.arange(240,len(t)))
@lru_cache(None)
def q2_stats(h,d,g):
 # Adapt Q2 fixed shift model to causal trailing flow, with common detector >=120s.
 a=D[h];i=a['idx'][::4];t=a['t'][i];u=a['O'][i]/3000
 q=a['q'][i-2*d]/5
 X=np.c_[np.ones(len(i)),t/600,u,np.maximum(u[:,None]-np.array([.15,.4,.7,.9]),0),q,q*q]
 y=a['y'][i,g]
 return X,y,X.T@X/len(i),X.T@y/len(i)
@lru_cache(None)
def lead(train,g):
 if args.lag_scheme=='q2_fixed':
  ds=range(6,61,6);loss=[]
  for d in ds:
   stats=[q2_stats(h,d,g) for h in train];A=sum(z[2] for z in stats);B=sum(z[3] for z in stats);pen=np.eye(A.shape[0])*.001*(len(train)-1);pen[0,0]=0
   loss.append(np.mean([np.abs(y-X@np.linalg.solve(A-aa+pen,B-bb)).mean() for X,y,aa,bb in stats]))
  return list(ds)[int(np.argmin(loss))]
 ds=np.arange(0,61,5);v=[]
 for h in train:
  a=D[h];i=a['idx'];y=a['y'][i,g]-a['y'][i-10,g];v.append([abs(np.corrcoef(a['br'][i-2*d]-a['br'][i-2*d-10],y)[0,1]) for d in ds])
 return int(ds[np.nanargmax(np.nanmean(v,0))])
@lru_cache(None)
def features(h,d):
 a=D[h];i=a['idx']-2*d;t=a['t'][i];O=a['O'][i];q=a['q'][i];u=np.maximum.accumulate(O+1.25*d*O/np.maximum(t,1));p=u/1000
 # process-only features, no spectrum
 X=np.c_[p,p*p,p*p*p,np.maximum(p[:,None]-np.arange(.5,4,.5),0)**2,q,a['q60'][i],O/np.maximum(t,1),t/100,q*p,a['q60'][i]*p]
 return u,X,a['S'][i],a['S32'][i]
@lru_cache(None)
def template(h,d,g):
 u=features(h,d)[0];y=uniform_filter1d(D[h]['y'][D[h]['idx'],g],24,mode='nearest');keep=np.r_[True,np.diff(u)>1e-7];return u[keep],y[keep]
G=np.linspace(0,5000,201)
@lru_cache(4096)
def prior(train,h,d,g):
 u=features(h,d)[0];V=np.array([np.interp(G,*template(j,d,g)) for j in train]);mu=V.mean(0);med=np.median(V,0);z=np.stack([D[j]['pref'] for j in train]);sd=np.maximum(z.std(0),1e-5);z=(z-z.mean(0))/sd;pt=np.clip((D[h]['pref']-np.stack([D[j]['pref'] for j in train]).mean(0))/sd,-3,3);s=max(mu.std(),1);coef=[]
 for v in V:
  A=np.c_[np.ones(len(mu)),mu/s];a,b=np.linalg.solve(A.T@A/len(mu)+np.eye(2)*.05,A.T@((v-mu)/s)/len(mu));coef.append([a,np.log(np.clip(1+b,.25,4))])
 a,b=Ridge(alpha=10).fit(z,coef).predict(pt[None])[0];sim=a*s+np.exp(np.clip(b,-.7,.7))*mu
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
 AA=np.array([x.T@x/len(x) for x in X]);bb=np.array([x.T@y/len(x) for x,y in zip(X,Y)]);reg=np.eye(xt.shape[1]);reg[0,0]=1e-8;outp={};a0=0 if base is None else prior(train,h,d,g)[base]
 for alpha in (.03,.3,3.):
  coef=np.linalg.solve(AA.mean(0)+alpha*reg,bb.mean(0));v=xt@coef
  for strength in ((1.,) if base is None else (.5,1.)):
   outp[f'{mode}_{base or "direct"}_ridge{alpha}_s{strength}']=a0+strength*v
  if base is None and mode in ('spectral','fusion'):
   # individually trained heat experts; global equal-weight mixture
   co=np.linalg.solve(AA+alpha*reg,bb[...,None]).squeeze(-1);outp[f'{mode}_experts_ridge{alpha}']=xt@co.mean(0)
 if base is None and mode in ('spectral','fusion'):
  xx=np.concatenate([x[:,1:] for x in X]);yy=np.concatenate(Y)
  for n in (3,8):
   mdl=PLSRegression(n_components=n,scale=False,max_iter=100).fit(xx,yy);outp[f'{mode}_PLS{n}']=mdl.predict(xt[:,1:]).ravel()
 return outp
@lru_cache(256)
def candidates(train,h,g):
 assert h not in train;d=lead(train,g);pr=prior(train,h,d,g);p={'process_template':pr['median'],'prior_sim':pr['sim'],'prior_fpcr':pr['fpcr']}
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
def feedback(v,y,d,interval,gain):
 # y indexed by detector clock, first scored detector = 120 s. At publish t, last observed detector timestamp <=t.
 k=np.arange(len(y));step=int(interval*2);last=((k-2*d)//step)*step;ok=last>=0;res=np.zeros(len(y));res[ok]=y[last[ok]]-v[last[ok]];return np.maximum(0,v+gain*res)

def main():
 metrics=[];inner=[];pred=[]
 for h in hs:
  tr=tuple(j for j in hs if j!=h)
  for g,gas in enumerate(('CO2','CO')):
   rows=[]
   for j in tr:
    sub=tuple(k for k in tr if k!=j);p,d=candidates(sub,j,g);y=D[j]['y'][D[j]['idx'],g]
    for n,v in p.items():
     f=family(n);rows.append(dict(model=n,family=f,inner_heat=letters[j],**score(y,v)))
     if f=='fusion':
      for gain in (.5,1.):
       nm=n+f'__fb{gain}';rows.append(dict(model=nm,family='feedback60',inner_heat=letters[j],**score(y,feedback(v,y,d,60,gain))))
   r=pd.DataFrame(rows);rank=r.groupby(['family','model']).r2.mean().reset_index();chosen=rank.sort_values('r2',ascending=False).drop_duplicates('family');inner.extend(r.assign(heat=letters[h],gas=gas).to_dict('records'))
   p,d=candidates(tr,h,g);y=D[h]['y'][D[h]['idx'],g];t=D[h]['t'][D[h]['idx']]
   for _,z in chosen.iterrows():
    n=z.model;f=z.family;v=p[n] if f!='feedback60' else feedback(p[n.split('__fb')[0]],y,d,60,float(n.split('__fb')[1]));metrics.append(dict(heat=letters[h],heat_id=h,gas=gas,family=f,model=n,inner_r2=z.r2,lead_s=d,n=len(y),**score(y,v)));pred.extend(zip([letters[h]]*len(y),[gas]*len(y),[f]*len(y),t-d,t,y,v))
  print(letters[h],flush=True);pd.DataFrame(metrics).to_csv(out/'metrics.csv',index=False)
 pd.DataFrame(inner).to_csv(out/'inner.csv.gz',index=False);pd.DataFrame(pred,columns=['heat','gas','family','publish_s','detector_s','observed','prediction']).to_csv(out/'predictions.csv.gz',index=False)
if __name__=='__main__':
 with threadpool_limits(limits=1):main()
