import numpy as np,pandas as pd,json
from pathlib import Path
from functools import lru_cache
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge,LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.ndimage import uniform_filter1d
from scipy.special import softmax
from threadpoolctl import threadpool_limits
import os
root=Path(os.environ['Q3_B_OUTPUT']); src=Path(os.environ['Q3_B_CACHE']); B='6B03206'; hs=tuple(p.stem for p in sorted(src.glob('6*.npz'))); tr=tuple(h for h in hs if h!=B); letters=dict(zip(hs,'ABCDEFGHIJ')); D={}
for h in hs:
 a=np.load(src/(h+'.npz')); O=a['O'].astype(float);O-=O[0];i=np.arange(len(O));t=i*.5
 def trailing(x,n):
  s=np.vstack([np.zeros((1,x.shape[1])),np.cumsum(x,axis=0)]);j=i+1;k=np.maximum(0,j-n);return (s[j]-s[k])/(j-k)[:,None]
 z=trailing(a['bands'].astype(float),30);c=trailing(a['contrast'].astype(float),30);br=z.mean(1);sh=z-br[:,None];k=np.maximum(0,i-20);q=(O-O[k])/np.maximum((i-k)*.5,.5); m=t<=60
 p=np.c_[br,sh.reshape(-1,8,4).mean(2),(sh-sh[m].mean(0)).reshape(-1,8,4).mean(2),c.reshape(-1,8,4).mean(2),q,O/np.maximum(t,1)]
 pref=np.r_[br[m].mean(),br[m].std(),q[m].mean(),q[m].std(),sh[m].mean(0).reshape(4,8).mean(1),c[m].mean()]
 D[h]=dict(O=O,t=t,p=p,pref=pref,br=br,idx=i[240:]);
 if h!=B:D[h]['y']=a['gas'].astype(float)
@lru_cache(None)
def lead(train,g):
 ds=np.arange(0,61,5);v=[]
 for h in train:
  a=D[h];i=a['idx'];y=a['y'][i,g]-a['y'][i-10,g];v.append([abs(np.corrcoef(a['br'][i-2*d]-a['br'][i-2*d-10],y)[0,1]) for d in ds])
 return int(ds[np.nanargmax(np.nanmean(v,0))])
@lru_cache(None)
def rep(h,d):
 a=D[h];j=a['idx']-2*d;O=a['O'][j];return np.maximum.accumulate(O+1.25*d*O/np.maximum(a['t'][j],1)),a['p'][j]
@lru_cache(None)
def template(h,d,g):
 u,p=rep(h,d);return u,p,uniform_filter1d(D[h]['y'][D[h]['idx'],g],24,mode='nearest')
def interp(u,v,x):
 k=np.r_[True,np.diff(u)>1e-7];u=u[k];v=v[k];return np.interp(x,u,v) if v.ndim==1 else np.stack([np.interp(x,u,v[:,j]) for j in range(v.shape[1])],1)
def pair(train,test,d,g):
 u,p=rep(test,d);sd=np.maximum(np.concatenate([rep(h,d)[1] for h in train]).std(0),1e-5);sp=np.maximum(np.stack([D[h]['pref'] for h in train]).std(0),1e-5);V=[];X=[]
 for h in train:
  uh,ph,yh=template(h,d,g);delta=np.minimum(((p-interp(uh,ph,u))/sd)**2,100);dp=np.minimum(((D[test]['pref']-D[h]['pref'])/sp)**2,100).mean();X.append(np.c_[delta[:,0],delta[:,1:9].mean(1),delta[:,9:17].mean(1),delta[:,17:25].mean(1),delta[:,25:].mean(1),np.full(len(u),dp),np.maximum(0,u-uh[-1])/1000,np.maximum(0,uh[0]-u)/1000,u/3000,(u/3000)**2]);V.append(interp(uh,yh,u))
 return np.array(V),np.array(X)
@lru_cache(None)
def router(train,d,g):
 X=[];Y=[]
 for h in train:
  sub=tuple(j for j in train if j!=h);v,x=pair(sub,h,d,g);u=rep(h,d)[0];y=D[h]['y'][D[h]['idx'],g];s=max(y.std(),1)
  for lo,hi in zip([0,600,1200,1800,2400,3000],[600,1200,1800,2400,3000,1e9]):
   m=(u>=lo)&(u<hi)
   if m.sum()>=8:X.extend(x[:,m].mean(1));Y.extend(np.log(.05+((v[:,m]-y[m])**2).mean(1)/s**2))
 sc=StandardScaler().fit(X);return sc,Ridge(alpha=30).fit(sc.transform(X),Y)
@lru_cache(None)
def base(train,test,g):
 d=lead(train,g);v,x=pair(train,test,d,g);med=np.median(v,0);out={'global':med};info={};grid=np.linspace(max(rep(h,d)[0][0] for h in train),min(rep(h,d)[0][-1] for h in train),60);sh=[]
 for h in train:
  u,_,y=template(h,d,g);a=interp(u,y,grid);sh.append((a-a.mean())/max(a.std(),.001))
 pref=np.stack([D[h]['pref'] for h in train]);sc=StandardScaler().fit(pref);pt=sc.transform(D[test]['pref'][None]);P=sc.transform(pref)
 for k in (2,3):
  lab=KMeans(n_clusters=k,n_init=10,random_state=123).fit_predict(sh);clf=LogisticRegression(C=.1,max_iter=1000).fit(P,lab);prob=clf.predict_proba(pt)[0];experts=[]
  for j in clf.classes_:
   m=lab==j;n=m.sum();w=n/(n+2);experts.append((1-w)*med+w*v[m].mean(0))
  out[f'type{k}']=prob@np.array(experts);info[str(k)]={'groups':{str(int(j)):[letters[h] for h,l in zip(train,lab) if l==j] for j in clf.classes_},'prob':prob.tolist()}
 sc,r=router(train,d,g);loss=np.clip(r.predict(sc.transform(x.reshape(-1,x.shape[-1]))).reshape(x.shape[:2]),-5,5)
 for T in (.5,1.,2.):
  w=softmax(-loss/T,axis=0);out[f'route{T}']=.5*med+.5*(w*v).sum(0)
  for k in (2,3):out[f'hybrid{k}_{T}']=.5*(out[f'type{k}']+out[f'route{T}'])
 info['route_weights_T1']={letters[h]:float(wt) for h,wt in zip(train,softmax(-loss,axis=0).mean(1))};return out,d,info
names=('global','type2','type3','route1.0','hybrid2_1.0','hybrid3_1.0')
@lru_cache(None)
def calibrate(train,test,g):
 vals={n:[] for n in names};P=[];scales=[]
 for h in train:
  sub=tuple(j for j in train if j!=h);pred=base(sub,h,g)[0];y=D[h]['y'][D[h]['idx'],g];s=max(y.std(),1);P.append(D[h]['pref']);scales.append(s)
  for n in names:
   v=pred[n];z=np.c_[np.ones(len(v)),v/s];delta=np.linalg.solve(z.T@z/len(y)+np.diag([.05,.05]),z.T@((y-v)/s)/len(y));vals[n].append([delta[0]*s,np.log(np.clip(1+delta[1],.25,4))])
 sc=StandardScaler().fit(P);X=sc.transform(P);pt=sc.transform(D[test]['pref'][None]);s=np.mean(scales);out={}
 for n in names:
  v=np.array(vals[n]);v[:,0]/=s;a,b=Ridge(alpha=20).fit(X,v).predict(pt)[0];out[n]=(float(a*s),float(np.exp(np.clip(b,-1.4,1.4))))
 return out

def candidates(train,test,g):
 p,d,info=base(train,test,g);p=dict(p);ab=calibrate(train,test,g)
 for n,(a,b) in ab.items():
  for w in (.5,1.):p[n+f'_ab{w}']=(1-w)*p[n]+w*(a+b*p[n])
 return {k:np.maximum(v,0) for k,v in p.items()},d,info,ab

def score(y,v):
 e=y-v;return dict(r2=1-(e@e)/np.sum((y-y.mean())**2),mae=np.abs(e).mean(),rmse=np.sqrt(np.mean(e**2)))
if __name__=='__main__':
 with threadpool_limits(limits=1):
  rows=[];vr=[];preds=[];meta={}
  for g,gas in enumerate(('CO2','CO')):
   for h in tr:
    sub=tuple(j for j in tr if j!=h);p,d,_,_=candidates(sub,h,g);y=D[h]['y'][D[h]['idx'],g]
    for n,v in p.items():vr.append(dict(gas=gas,heat=letters[h],model=n,**score(y,v)))
    print('inner',gas,letters[h],d,flush=True)
   df=pd.DataFrame(vr);rank=df[df.gas==gas].groupby('model')[['r2','mae','rmse']].mean().sort_values('r2',ascending=False);rank.to_csv(root/f'{gas}_inner_summary.csv');chosen=rank.index[0];p,d,info,ab=candidates(tr,B,g);meta[gas]=dict(lead=d,selected=chosen,info=info,alpha_beta=ab)
   # B outcome not loaded before prediction and selection.
   y=np.load(src/(B+'.npz'))['gas'][D[B]['idx'],g].astype(float);det=D[B]['t'][D[B]['idx']]
   for n,v in p.items():
    rows.append(dict(gas=gas,model=n,selected=n==chosen,lead_s=d,n=len(y),**score(y,v)));preds.extend(zip([gas]*len(y),[n]*len(y),det-d,det,y,v))
   pd.DataFrame(rows).to_csv(root/'B_metrics.csv',index=False);pd.DataFrame(vr).to_csv(root/'inner_folds.csv',index=False);pd.DataFrame(preds,columns=['gas','model','publish_s','detector_s','observed','prediction']).to_csv(root/'B_predictions.csv',index=False);(root/'metadata.json').write_text(json.dumps(meta,indent=2));print('SELECTED',gas,chosen,score(y,p[chosen]),flush=True)
