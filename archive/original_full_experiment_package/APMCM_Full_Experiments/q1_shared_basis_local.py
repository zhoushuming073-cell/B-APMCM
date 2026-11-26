import os,glob,re,json,time
import numpy as np,pandas as pd
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import median_filter
from sklearn.decomposition import PCA
from sklearn.utils.extmath import randomized_svd
R='/mnt/data/APMCM_Full_Experiments'; files=sorted(glob.glob('/mnt/data/heats_*.npz'))
D={}
for f in files:
 z=np.load(f);h=os.path.basename(f)[6:-4];D[h]={k:z[k] for k in z.files}
# common grid: exact 6B grid restricted to overlap; C spectra are interpolated by wavelength, never by column index
wB=D['6B03205']['wl']; wC=D['6C02057']['wl']; grid=wB[(wB>=wC.min())&(wB<=wC.max())]
print('common grid',len(grid),grid[0],grid[-1],flush=True)
def harmonize(h):
 z=D[h];X=z['X'].astype(np.float32);w=z['wl']
 if np.array_equal(w,grid): return X.copy()
 # Interpolate only within measured support; saturation mask interpolated conservatively.
 out=np.empty((len(X),len(grid)),np.float32)
 for i,row in enumerate(X): out[i]=np.interp(grid,w,row).astype(np.float32)
 return out
H={h:harmonize(h) for h in D}
# save harmonized spectra and metadata
np.savez_compressed(R+'/harmonized_spectra.npz',grid=grid,**{h:H[h] for h in H})
# candidate peak counts in selected rows, excluding saturated neighborhoods
peak_counts=[]
for h,X in H.items():
 z=D[h];sat=(z['X']>=65535).any(axis=1)
 cnt=0; rows=0
 for i in range(20,len(X)-20,20):
  if sat[i]: continue
  p,props=find_peaks(X[i],prominence=max(100,np.nanmedian(X[i])*.08),width=(1,20),distance=4)
  p=p[(p>50)&(p<len(grid)-50)]
  cnt+=len(p);rows+=1
 peak_counts.append({'heat':h,'rows_sampled':rows,'peak_candidates':cnt,'per_row':cnt/max(rows,1)})
print('peak counts',pd.DataFrame(peak_counts).to_string(index=False),flush=True)
# training samples for PCA, exclude saturated rows and evenly sample max 80 per heat
train_cache={}
for h,X in H.items():
 sat=(D[h]['X']>=65535).any(axis=1)
 ids=np.linspace(0,len(X)-1,min(80,len(X))).astype(int);ids=ids[~sat[ids]]
 if len(ids)<20: ids=np.where(~sat)[0][::max(1,np.sum(~sat)//80)][:80]
 rows=X[ids].astype(np.float32)
 log=np.log1p(np.maximum(rows,0)); base=np.median(log,axis=1,keepdims=True); shape=log-base
 train_cache[h]=shape
# cross-validation masks: test 30 rows per heat, peak-centred widths 3/7/11, fixed RNG
rng=np.random.default_rng(20261126)
rows=[]; details=[]
# helpers
def fit_basis(train_heats, rank):
 A=np.vstack([train_cache[h] for h in train_heats])
 mu=A.mean(axis=0)
 # randomized SVD
 U,s,Vt=randomized_svd(A-mu,n_components=rank,n_iter=5,random_state=13)
 return mu,Vt

def pca_predict(row,mask,mu,Vt):
 y=np.log1p(np.maximum(row,0)); lev=np.median(y[~mask]); shape=y-lev
 obs=~mask
 B=Vt[:,obs].T
 coef=np.linalg.lstsq(B,shape[obs]-mu[obs],rcond=None)[0]
 predshape=mu+coef@Vt
 return np.expm1(np.maximum(0,predshape+lev))

def linear_predict(row,mask):
 idx=np.arange(len(row));obs=~mask
 return np.interp(idx,idx[obs],row[obs])

def temporal_predict(X,i,mask,j,width):
 z=np.log1p(np.maximum(X,0));half=width//2;left=max(0,j-half-12);right=min(X.shape[1],j+half+13)
 sh=np.r_[np.arange(left,max(left,j-half)),np.arange(min(right,j+half+1),right)]
 sh=sh[~mask[sh]]
 vals=[];weights=[]
 for k in [i-1,i+1]:
  if 0<=k<len(X) and len(sh)>=6:
   gain=np.median(z[i,sh]-z[k,sh]);vals.append(z[k]+gain);weights.append(1.0)
 if not vals:return np.expm1(z[i])
 return np.expm1(np.average(np.stack(vals),axis=0,weights=weights))

# For each held-out heat, train pooled basis excluding it; also same-family pooled basis, and within-heat early basis.
for h,X in H.items():
 z=D[h]; satrow=(z['X']>=65535).any(axis=1)
 candidate_rows=[]
 for i in range(25,len(X)-25,25):
  if satrow[i]: continue
  p,props=find_peaks(X[i],prominence=max(100,np.median(X[i])*.08),width=(1,20),distance=4)
  p=[int(v) for v in p if 60<=v<len(grid)-60]
  if p:
   rng.shuffle(p);candidate_rows.extend([(i,j) for j in p[:2]])
 rng.shuffle(candidate_rows);candidate_rows=candidate_rows[:60]
 if len(candidate_rows)<15: print('low peak candidates',h,len(candidate_rows),flush=True)
 all_train=[g for g in H if g!=h]
 same_train=[g for g in all_train if g[:2]==h[:2]]
 # For within-heat reference, train basis on first 55% rows, excluding saturation
 ids=np.arange(int(.55*len(X)));ids=ids[~satrow[ids]][::max(1,np.sum(~satrow[ids])//80)][:80]
 own=np.log1p(X[ids])-np.median(np.log1p(X[ids]),axis=1,keepdims=True)
 # Build all ranks and training sets
 basis={}
 for label,hs in [('pooled9',all_train),('same_family',same_train)]:
  if not hs: continue
  for rank in [3,6,12]: basis[(label,rank)]=fit_basis(hs,rank)
 for rank in [3,6,12]:
  mu=own.mean(axis=0);_,ss,Vt=randomized_svd(own-mu,n_components=rank,n_iter=4,random_state=13);basis[('own_early',rank)]=(mu,Vt)
 # each candidate and width
 for i,j in candidate_rows:
  for width in [3,7,11]:
   mask=np.zeros(len(grid),bool);half=width//2;mask[j-half:j+half+1]=True
   # don't evaluate saturated or edge
   if np.any(X[i,mask]>=65535): continue
   truth=X[i,mask]
   preds={'linear':linear_predict(X[i],mask)[mask]}
   tp=temporal_predict(X,i,mask,j,width);preds['temporal_local_shoulder'] = tp[mask]
   for key,(mu,Vt) in basis.items():
    pred=pca_predict(X[i],mask,mu,Vt);preds[f'{key[0]}_r{key[1]}']=pred[mask]
   # metrics; peak height and local excess area against local shoulder baseline
   shoulder=np.r_[max(0,j-half-8):j-half, j+half+1:min(len(grid),j+half+9)]
   shoulder=shoulder[~mask[shoulder]]
   baseline=np.median(X[i,shoulder]) if len(shoulder) else np.median(truth)
   true_area=np.sum(np.maximum(truth-baseline,0));true_height=np.max(truth)-baseline
   for method,pred in preds.items():
    mae=float(np.mean(np.abs(pred-truth)));rmse=float(np.sqrt(np.mean((pred-truth)**2)))
    pa=np.sum(np.maximum(pred-baseline,0));ph=np.max(pred)-baseline
    details.append({'heat':h,'row':i+1,'peak_idx':j,'width':width,'method':method,'mae':mae,'rmse':rmse,'peak_height_rel_err':abs(ph-true_height)/max(abs(true_height),1),'peak_area_rel_err':abs(pa-true_area)/max(abs(true_area),1),'true_peak_height':true_height,'true_peak_area':true_area})
 print('finished LOHO',h,'peaks',len(candidate_rows),flush=True)
res=pd.DataFrame(details);res.to_csv(R+'/q1_peak_masking_local_long.csv',index=False)
summary=res.groupby(['method','width']).agg(n=('mae','size'),MAE=('mae','mean'),RMSE=('rmse','mean'),height_rel=('peak_height_rel_err','median'),area_rel=('peak_area_rel_err','median')).reset_index()
summary.to_csv(R+'/q1_peak_masking_local_summary.csv',index=False)
print(summary.pivot(index='method',columns='width',values='MAE').round(1).to_string(),flush=True)
print('area errors',summary.pivot(index='method',columns='width',values='area_rel').round(3).to_string(),flush=True)
# Cross-family calibration impact: compare pooled all9 vs same-family vs own basis per heat and family
