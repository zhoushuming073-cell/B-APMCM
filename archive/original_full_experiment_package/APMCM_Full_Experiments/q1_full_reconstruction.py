from pathlib import Path
import os,glob,sys
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d as uf,maximum_filter,median_filter
R=Path('/mnt/data/APMCM_Full_Experiments');outdir=R/'cleaned_scenarios';outdir.mkdir(exist_ok=True)

def fit_scales(train):
 sc={}
 for dt in [1,3,5]:
  za=uf(np.log1p(train),dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5;rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2;r=(rp+rn)/2
   sc[(dt,dw)]=np.maximum(1.4826*np.median(abs(r[15:-15]),axis=0),.002)
 return sc

def score(a,sc):
 z=np.log1p(a);multi=np.zeros(a.shape,np.float32)
 for dt in [1,3,5]:
  za=uf(z,dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5;rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2;r=(rp+rn)/2;ag=np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0);v=.65*abs(r)+.35*ag;base=sc[(dt,dw)];w=uf(abs(r),41,axis=1);off2=dw//2+21;local=(np.roll(w,off2,axis=1)+np.roll(w,-off2,axis=1))/2;sl=v/np.maximum(base,1.4826*local);np.maximum(multi,maximum_filter(sl,size=(dt,dw)),out=multi)
 return multi

def options(a,mask):
 z=np.log1p(a);p=np.full(a.shape,np.nan,np.float32);q=p.copy()
 for i in np.where(mask.any(axis=1))[0]:
  js=np.flatnonzero(mask[i]);runs=np.split(js,np.where(np.diff(js)>1)[0]+1)
  for g in runs:
   if len(g)>25:continue
   l,r=g[0],g[-1];sh=np.r_[np.arange(max(0,l-15),l),np.arange(r+1,min(a.shape[1],r+16))]
   for direction,out in [(-1,p),(1,q)]:
    for d in range(1,13):
     k=i+direction*d
     if not 0<=k<len(a):break
     if mask[k,g].any():continue
     ids=sh[~(mask[i,sh]|mask[k,sh])]
     if len(ids)<8:continue
     out[i,g]=z[k,g]+np.median(z[i,ids]-z[k,ids]);break
 v=(p+q)/2;diff=np.minimum(abs(z-p),abs(z-q));agree=(z-p)*(z-q)>0;gap=abs(p-q)
 return v,diff,agree,gap

def repair_sat(a,sat):
 out=a.copy();done=np.zeros(a.shape,bool);n,m=a.shape
 for i in np.where(sat.any(axis=1))[0]:
  js=np.flatnonzero(sat[i]);runs=np.split(js,np.where(np.diff(js)>1)[0]+1)
  for g in runs:
   l,r=int(g[0]),int(g[-1]);should=np.r_[np.arange(max(0,l-15),l),np.arange(r+1,min(m,r+16))]
   if len(should)<8:continue
   vals=[]
   for k in [i-1,i+1]:
    if not 0<=k<n or sat[k,g].any():continue
    ids=should[~(sat[i,should]|sat[k,should])]
    if len(ids)<8:continue
    gain=np.median(np.log1p(a[i,ids])-np.log1p(a[k,ids]));vals.append(np.log1p(a[k,g])+gain)
   if vals:
    out[i,g]=np.expm1(np.mean(vals,axis=0));done[i,g]=True
 return out,done

rows=[]
files=sorted(glob.glob('/mnt/data/heats_*.npz'))
if len(sys.argv)>1: files=[f for f in files if os.path.basename(f)[6:-4] in sys.argv[1:]]
for f in files:
 z=np.load(f);h=os.path.basename(f)[6:-4];X=z['X'].astype(np.float32);n,m=X.shape;sat=X>=65535
 cut=int(.55*n);train=X[:cut];sc=fit_scales(train);sc0=score(train,sc);valid=np.ones(train.shape,bool);valid[:15]=valid[-15:]=False;valid[:,:40]=valid[:,-40:]=False;th=float(np.quantile(sc0[valid],.999))
 fullscore=score(X,sc);vm=np.ones(X.shape,bool);vm[:15]=vm[-15:]=False;vm[:,:40]=vm[:,-40:]=False
 screen=(fullscore>th)&vm&~sat
 v,diff,agree,gap=options(X,screen);approved=screen&agree&(diff>.12)&(diff>2*gap)&np.isfinite(v)
 scen=X.copy();scen[approved]=np.expm1(v[approved])
 sat_scen,sat_done=repair_sat(scen,sat);scen=sat_scen
 np.savez_compressed(outdir/f'{h}_scenario.npz',raw=X,scenario=scen,wl=z['wl'],seq=z['seq'],score=fullscore,screen=screen,approved=approved,saturated=sat,saturation_repaired=sat_done)
 ii,jj=np.where(screen);ai,aj=np.where(approved);si,sj=np.where(sat);di,dj=np.where(sat_done)
 rows.append({'heat':h,'threshold':th,'screen_candidates':len(ii),'screen_rate':len(ii)/vm.sum(),'approved_anomaly_repairs':len(ai),'approved_rate':len(ai)/vm.sum(),'saturation_points':len(si),'saturation_repaired':len(di),'saturation_unresolved':len(si)-len(di),'changed_total':int(approved.sum()+sat_done.sum()),'large_runs_skipped':int(np.sum(screen&~approved))})
 # candidate detail top 500 by score
 if len(ii):
  ix=np.argsort(fullscore[ii,jj])[::-1][:500];pd.DataFrame({'sequence':ii[ix]+1,'time_s':(ii[ix])*.5,'wavelength_nm':z['wl'][jj[ix]],'raw':X[ii[ix],jj[ix]],'score':fullscore[ii[ix],jj[ix]],'approved_repair':approved[ii[ix],jj[ix]],'scenario':scen[ii[ix],jj[ix]],'saturated':sat[ii[ix],jj[ix]]}).to_csv(R/f'{h}_candidate_top500.csv',index=False)
 print('scenario',h,'screen',len(ii),'approved',len(ai),'sat',len(si),'sat repaired',len(di),flush=True)
pd.DataFrame(rows).to_csv(R / ('q1_full_reconstruction_summary_'+('_'.join(sys.argv[1:]) if len(sys.argv)>1 else 'all')+'.csv'),index=False)
print(pd.DataFrame(rows).to_string(index=False))
