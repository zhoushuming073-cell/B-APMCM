from pathlib import Path
import os,sys
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d as uf, maximum_filter,median_filter
R=Path('/mnt/data/APMCM_Full_Experiments');h=sys.argv[1];z=np.load('/mnt/data/heats_'+h+'.npz');X=z['X'].astype(np.float32);n,m=X.shape
shapes=[(1,1),(1,3),(1,5),(3,1),(3,3),(5,5)]
def reg(e):i,j,dt,dw,*_=e;return slice(i-dt//2,i+dt//2+1),slice(j-dw//2,j+dw//2+1)
def inject(a,seed,per=4):
 rg=np.random.default_rng(seed);b=a.copy();y=np.zeros(a.shape,bool);ev=[]
 for dt,dw in shapes:
  for s in [1,-1]:
   for amp in [.05,.15,.4]:
    for _ in range(per):
     for tr in range(3000):
      i=int(rg.integers(20,len(a)-20));j=int(rg.integers(45,m-45));e=(i,j,dt,dw,s,amp);ss=reg(e)
      if y[i-16:i+17,j-35:j+36].any() or a[ss].min()<300 or a[ss].max()>45000:continue
      break
     else:continue
     b[ss]=a[ss]*(1+s*amp);y[ss]=True;ev.append(e)
 return b,y,ev
def valid(a):
 v=np.ones(a.shape,bool);v[:15]=v[-15:]=False;v[:,:40]=v[:,-40:]=False;return v
def scales_fit(train):
 sc={}
 for dt in [1,3,5]:
  za=uf(np.log1p(train),dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5;rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2;r=(rp+rn)/2
   sc[(dt,dw)]=np.maximum(1.4826*np.median(abs(r[15:-15]),axis=0),.002)
 z=np.log1p(train);p=z-np.roll(z,1,axis=0);q=z-np.roll(z,-1,axis=0);r=((p-median_filter(p,size=(1,21)))+(q-median_filter(q,size=(1,21))))/2;sc['prior']=np.maximum(1.4826*np.median(abs(r[4:-4]),axis=0),.004)
 return sc
def score(a,sc):
 z=np.log1p(a);multi=np.zeros(a.shape,np.float32)
 for dt in [1,3,5]:
  za=uf(z,dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5;rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2;r=(rp+rn)/2;ag=np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0);v=.65*abs(r)+.35*ag;base=sc[(dt,dw)];w=uf(abs(r),41,axis=1);off2=dw//2+21;local=(np.roll(w,off2,axis=1)+np.roll(w,-off2,axis=1))/2;sl=v/np.maximum(base,1.4826*local);np.maximum(multi,maximum_filter(sl,size=(dt,dw)),out=multi)
 return multi
# two independent local shoulder predictions for candidate regions
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
 return z,v,diff,agree,gap
cut=int(.55*n);train=X[:cut];test0=X[min(cut+40,n-30):]
if len(test0)>350:test0=test0[-350:]
if len(test0)<120:test0=X[-350:]
sc=scales_fit(train);co=score(train,sc);cv=valid(train);th=float(np.quantile(co[cv],.999))
ci,cy,ce=inject(train,20261126+sum(map(ord,h)));ti,ty,te=inject(test0,20261226+sum(map(ord,h)))
cm=(score(ci,sc)>th)&cv;tm=(score(ti,sc)>th)&valid(test0)
copt=options(ci,cm);topt=options(ti,tm)
rows=[]
for af in [.025,.05,.08,.12,.18]:
 for gf in [.5,1,2]:
  for part,raw,true,mask,opt,vm in [('cal',ci,cy,cm,copt,cv),('test',ti,ty,tm,topt,valid(test0))]:
   z0,v,d,ag,gap=opt;ok=mask&ag&(d>af)&(d>gf*gap)&np.isfinite(v);pred=raw.copy();pred[ok]=np.expm1(v[ok]);err=abs(pred-(train if part=='cal' else test0));e0=abs(raw-(train if part=='cal' else test0));rows.append({'heat':h,'partition':part,'absfloor':af,'gapfactor':gf,'screen_point_recall':float(mask[true].mean()),'repair_point_recall':float(ok[true].mean()),'event_recall':float(np.mean([ok[reg(e)].any() for e in (ce if part=='cal' else te)])),'injected_mae_before':float(e0[true].mean()),'injected_mae_after':float(err[true].mean()),'collateral_mae':float(err[~true & vm].mean()),'total_mae_before':float(e0[vm].mean()),'total_mae_after':float(err[vm].mean()),'changed_uninjected':int((ok&~true).sum())})
res=pd.DataFrame(rows);res.to_csv(R/f'gate_grid_{h}.csv',index=False)
cal=res[res.partition=='cal'];best=cal.sort_values('total_mae_after').iloc[0];sel=res[(res.absfloor==best.absfloor)&(res.gapfactor==best.gapfactor)&(res.partition=='test')].iloc[0]
print('HEAT',h,'threshold',th,'best gate',best.absfloor,best.gapfactor,'cal',best.total_mae_before,best.total_mae_after,'test',sel.total_mae_before,sel.total_mae_after,'test injected',sel.injected_mae_before,sel.injected_mae_after,'recall',sel.repair_point_recall,'event',sel.event_recall,flush=True)
