import os,glob
import numpy as np,pandas as pd
from scipy.signal import find_peaks
from sklearn.utils.extmath import randomized_svd
R='/mnt/data/APMCM_Full_Experiments';D={}
for f in glob.glob('/mnt/data/heats_*.npz'):
 z=np.load(f);h=os.path.basename(f)[6:-4];D[h]={k:z[k] for k in z.files}
train_cache={}
for h,z in D.items():
 X=z['X'].astype(np.float32);sat=(X>=65535).any(1);ids=np.linspace(0,len(X)-1,min(80,len(X))).astype(int);ids=ids[~sat[ids]];A=np.log1p(X[ids]);train_cache[h]=A-np.median(A,axis=1,keepdims=True)
def fit(hs,rank):
 A=np.vstack([train_cache[h] for h in hs]);mu=A.mean(0);_,s,Vt=randomized_svd(A-mu,n_components=rank,n_iter=5,random_state=13);return mu,Vt
def pred(row,mask,mu,Vt):
 y=np.log1p(row);lev=np.median(y[~mask]);sh=y-lev;B=Vt[:,~mask].T;coef=np.linalg.lstsq(B,sh[~mask]-mu[~mask],rcond=None)[0];return np.expm1(np.maximum(0,mu+coef@Vt+lev))
def localtemp(X,i,j,width,mask):
 z=np.log1p(X);half=width//2;sh=np.r_[np.arange(max(0,j-half-12),j-half),np.arange(j+half+1,min(X.shape[1],j+half+13))];sh=sh[~mask[sh]];v=[]
 for k in [i-1,i+1]:
  if 0<=k<len(X) and len(sh)>=6:
   gain=np.median(z[i,sh]-z[k,sh]);v.append(z[k]+gain)
 return np.expm1(np.mean(v,axis=0)) if v else X[i]
rows=[]
for h,z in D.items():
 X=z['X'].astype(np.float32);wl=z['wl'];n,m=X.shape;sat=(X>=65535).any(1);family=h[:2];hs=[g for g in D if g!=h and g[:2]==family]
 ids=np.arange(int(.55*n));ids=ids[~sat[ids]][::max(1,np.sum(~sat[ids])//80)][:80];own=np.log1p(X[ids]);own=own-np.median(own,axis=1,keepdims=True);mu0=own.mean(0);_,_,V0=randomized_svd(own-mu0,n_components=12,n_iter=4,random_state=13)
 basis={r:fit(hs,r) for r in [6,12,18]}
 cand=[];rng=np.random.default_rng(20261127+sum(map(ord,h)))
 for i in range(25,n-25,25):
  if sat[i]:continue
  p,pr=find_peaks(X[i],prominence=max(100,np.median(X[i])*.08),width=(1,20),distance=4);p=[int(j) for j in p if 60<=j<m-60]
  rng.shuffle(p);cand.extend([(i,j) for j in p[:2]])
 rng.shuffle(cand);cand=cand[:60]
 for i,j in cand:
  for width in [3,7,11]:
   mask=np.zeros(m,bool);half=width//2;mask[j-half:j+half+1]=True
   if np.any(X[i,mask]>=65535):continue
   truth=X[i,mask];should=np.r_[np.arange(j-half-8,j-half),np.arange(j+half+1,j+half+9)];base=np.median(X[i,should]);ta=np.sum(np.maximum(truth-base,0))
   preds={'linear':np.interp(np.arange(m),np.flatnonzero(~mask),X[i,~mask])[mask],'temporal_local':localtemp(X,i,j,width,mask)[mask]}
   for r,(mu,Vt) in basis.items():preds[f'same_native_r{r}']=pred(X[i],mask,mu,Vt)[mask]
   preds['own_native_r12']=pred(X[i],mask,mu0,V0)[mask]
   for method,pr in preds.items():
    ph=np.max(pr)-base;pa=np.sum(np.maximum(pr-base,0));rows.append({'heat':h,'family':family,'width':width,'method':method,'MAE':np.mean(abs(pr-truth)),'peak_height_rel_error':abs(ph-(np.max(truth)-base))/max(abs(np.max(truth)-base),1),'peak_area_rel_error':abs(pa-ta)/max(abs(ta),1)})
 print('native done',h,'peaks',len(cand),flush=True)
df=pd.DataFrame(rows);df.to_csv(R+'/q1_native_basis_long.csv',index=False);s=df.groupby(['method','width']).agg(n=('MAE','size'),MAE=('MAE','mean'),height_rel=('peak_height_rel_error','median'),area_rel=('peak_area_rel_error','median')).reset_index();s.to_csv(R+'/q1_native_basis_summary.csv',index=False);print(s.pivot(index='method',columns='width',values='MAE').round(1).to_string());print(s.pivot(index='method',columns='width',values='area_rel').round(3).to_string())
