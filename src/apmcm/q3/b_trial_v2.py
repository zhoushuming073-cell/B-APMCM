"""B-only stage-dependent MoE and affine calibration; B outcome loaded after selection."""
import os,json
from pathlib import Path
from functools import lru_cache
import numpy as np,pandas as pd
from scipy.special import softmax
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from threadpoolctl import threadpool_limits
from . import b_trial_v1 as old
ROOT=Path(os.environ['Q3_B_OUTPUT'])
D,B,tr,letters=old.D,old.B,old.tr,old.letters
@lru_cache(None)
def pack(train,h,g):
 d=old.lead(train,g);u,p=old.rep(h,d);v,x=old.pair(train,h,d,g)
 pref=np.stack([D[j]['pref'] for j in train]);mu=pref.mean(0);sd=np.maximum(pref.std(0),1e-5)
 P=np.clip((pref-mu)/sd,-3,3);ph=np.clip((D[h]['pref']-mu)/sd,-3,3)
 grid=np.linspace(max(old.rep(j,d)[0][0] for j in train),min(old.rep(j,d)[0][-1] for j in train),60)
 shapes=[];pp=[]
 for j in train:
  uj,pj,yj=old.template(j,d,g);yy=old.interp(uj,yj,grid);shapes.append((yy-yy.mean())/max(yy.std(),.01));pp.append(old.interp(uj,pj,u))
 med=np.median(v,0);out={'global':med};info={}
 dist=np.minimum(x[:,:,:5],16).mean(2);sm=softmax(-dist/2,axis=0)
 out['local']=.5*med+.5*(sm*v).sum(0)
 for K in (2,3):
  lab=KMeans(n_clusters=K,n_init=10,random_state=123).fit_predict(shapes)
  clf=LogisticRegression(C=.1,max_iter=1000).fit(P,lab);probs=clf.predict_proba(ph[None])[0];experts=[];dyn=[]
  for k in clf.classes_:
   ids=lab==k;n=ids.sum();shrink=n/(n+2);experts.append((1-shrink)*med+shrink*v[ids].mean(0));dyn.append(sm[ids].sum(0))
  ex=np.array(experts);dyn=np.array(dyn);probs2=probs[:,None]*np.sqrt(np.maximum(dyn,1e-8));probs2/=probs2.sum(0)
  out[f'type{K}']=probs@ex;out[f'dyn{K}']=(probs2*ex).sum(0)
  info[str(K)]={'groups':{str(int(k)):[letters[j] for j,l in zip(train,lab) if l==k] for k in clf.classes_},'p_prefix':probs.tolist(),'p_dynamic_mean':probs2.mean(1).tolist()}
 rel=p-np.median(pp,axis=0)
 zrel=np.c_[rel[:,0],rel[:,1:9].reshape(-1,4,2).mean(2),rel[:,9:17].reshape(-1,4,2).mean(2),rel[:,17:25].mean(1),rel[:,-2:]]
 return out,u,np.c_[zrel,p[:,-2:],p[:,0]],info,d

def phase(u):
 z=u/1000;centers=np.array([.5,1.5,2.5,3.5]);b=np.maximum(1-np.abs(z[:,None]-centers),0);i=b.sum(1)==0;b[i,np.argmin(np.abs(z[i,None]-centers),axis=1)]=1
 return b/b.sum(1,keepdims=True)
@lru_cache(None)
def inputs(train,g):
 rows=[]
 for j in train:
  sub=tuple(k for k in train if k!=j);p,u,x,_,_=pack(sub,j,g);y=D[j]['y'][D[j]['idx'],g];rows.append((j,p,u,x,y))
 return rows

def candidates(train,h,g):
 bases,u,x,info,d=pack(train,h,g);pred=dict(bases);rows=inputs(train,g);scale=max(np.mean([z[-1].std() for z in rows]),1.)
 pv=np.stack([D[j]['pref'] for j,_,_,_,_ in rows]);pm=pv.mean(0);ps=np.maximum(pv.std(0),1e-5);P=np.clip((pv-pm)/ps,-2.5,2.5);pt=np.clip((D[h]['pref']-pm)/ps,-2.5,2.5)
 dist=np.mean((P-pt)**2,axis=1);w=softmax(-dist/1.5);trust=float(np.exp(-max(0,dist.min()-1)/2))
 raw=[z[3] for z in rows];fm=np.mean([z.mean(0) for z in raw],axis=0);fs=np.maximum(np.sqrt(np.mean([((z-fm)**2).mean(0) for z in raw],axis=0)),1e-5)
 for name in ('global','type3','dyn3','local'):
  v=bases[name];center=float(np.mean([z[1][name].mean() for z in rows]));co=[];cp=[]
  for j,basesj,uj,xx,y in rows:
   b=basesj[name];z=(b-center)/scale;e=(y-b)/scale;H=phase(uj);C=np.c_[np.ones(len(b)),z];A=np.c_[H,H*z[:,None]]
   def est(M,reg):return np.linalg.solve(M.T@M/len(e)+np.eye(M.shape[1])*reg,M.T@e/len(e))
   co.append(est(C,.2));cp.append(est(A,.15))
  V=(v-center)/scale;H=phase(u);C=np.c_[np.ones(len(v)),V];A=np.c_[H,H*V[:,None]]
  for label,coeff,mat in [('static',co,C),('phase',cp,A)]:
   coeff=np.array(coeff);c=(.5*coeff.mean(0)+.5*(w@coeff))*trust
   for mix in (.25,.5,1.):pred[f'{name}_{label}_{mix}']=v+mix*scale*(mat@c)
  def design(b,uu,xx):
   q=np.clip((xx-fm)/fs,-3,3);H=phase(uu);bv=np.clip((b-center)/scale,-4,4);return np.c_[H,H*bv[:,None],q,q*(uu/3000-.5)[:,None]]
  X=[design(b[name],uj,xx) for _,b,uj,xx,y in rows];Y=[(y-b[name])/scale for _,b,uj,xx,y in rows];xt=design(v,u,x)
  AA=np.mean([z.T@z/len(z) for z in X],axis=0);bb=np.mean([z.T@y/len(y) for z,y in zip(X,Y)],axis=0)
  for reg in (.3,1.,3.):
   a=np.linalg.solve(AA+reg*np.eye(len(bb)),bb)
   for mix in (.5,1.):pred[f'{name}_resid_{reg}_{mix}']=v+mix*scale*np.clip(xt@a,-1.5,1.5)
 return {k:np.maximum(v,0) for k,v in pred.items()},d,info

def main():
 rows=[];predictions=[];results=[];locked={}
 for g,gas in enumerate(('CO2','CO')):
  for h in tr:
   sub=tuple(j for j in tr if j!=h);p,d,_=candidates(sub,h,g);y=D[h]['y'][D[h]['idx'],g]
   rows.extend(dict(gas=gas,heat=letters[h],model=n,**old.score(y,v)) for n,v in p.items());print('inner',gas,letters[h],flush=True)
  f=pd.DataFrame(rows);rank=f[f.gas==gas].groupby('model')[['r2','mae','rmse']].mean().sort_values('r2',ascending=False);rank.to_csv(ROOT/f'{gas}_inner.csv');chosen=rank.index[0];p,d,info=candidates(tr,B,g);locked[gas]=(p,d,chosen,info);print('LOCKED',gas,chosen,rank.iloc[0].to_dict(),flush=True)
  pd.DataFrame(rows).to_csv(ROOT/'inner_folds.csv',index=False)
 (ROOT/'locked.json').write_text(json.dumps({gas:{'chosen':a[2],'lead':a[1],'groups':a[3]} for gas,a in locked.items()},indent=2))
 for g,gas in enumerate(('CO2','CO')):
  p,d,chosen,info=locked[gas];y=np.load(old.src/(B+'.npz'))['gas'][D[B]['idx'],g];det=D[B]['t'][D[B]['idx']]
  for n,v in p.items():
   results.append(dict(gas=gas,model=n,selected=n==chosen,lead_s=d,n=len(y),**old.score(y,v)));predictions.extend(zip([gas]*len(y),[n]*len(y),det-d,det,y,v))
 pd.DataFrame(results).to_csv(ROOT/'B_metrics.csv',index=False);pd.DataFrame(predictions,columns=['gas','model','publish_s','detector_s','observed','prediction']).to_csv(ROOT/'B_predictions.csv.gz',index=False,compression={'method':'gzip','mtime':0});print(pd.DataFrame(results).query('selected').to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=1):main()
