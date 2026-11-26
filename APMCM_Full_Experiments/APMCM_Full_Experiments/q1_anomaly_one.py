from pathlib import Path
import os,sys
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d as uf, maximum_filter,median_filter
from scipy.signal import find_peaks
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
def template_scores(a,scales):
 z=np.log1p(a);multi=np.zeros(a.shape,np.float32)
 for dt in [1,3,5]:
  za=uf(z,dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5
   rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2
   r=(rp+rn)/2;ag=np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0);v=.65*abs(r)+.35*ag
   base=scales[(dt,dw)];w=uf(abs(r),41,axis=1);off2=dw//2+21;local=(np.roll(w,off2,axis=1)+np.roll(w,-off2,axis=1))/2
   sl=v/np.maximum(base,1.4826*local);np.maximum(multi,maximum_filter(sl,size=(dt,dw)),out=multi)
 p=z-np.roll(z,1,axis=0);q=z-np.roll(z,-1,axis=0);rp=p-median_filter(p,size=(1,21));rn=q-median_filter(q,size=(1,21));r=(rp+rn)/2;ag=np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0)
 prior=(.65*abs(r)+.35*ag)/np.maximum(np.maximum(1.4826*median_filter(abs(r),size=(1,21)),scales['prior']),.004)
 hamp=abs(z-median_filter(z,size=(9,1)))/scales['hampel']
 return {'multiscale_local':multi,'prior_soft':prior,'hampel':hamp}
# 55% training; late test slice (max 350 frames) to reduce repeated temporal overlap
cut=int(.55*n);train=X[:cut];test0=X[min(cut+40,n-30):]
if len(test0)>350: test0=test0[-350:]
if len(test0)<120: test0=X[-350:]
# scale calibration
scales={}
for dt in [1,3,5]:
 za=uf(np.log1p(train),dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
 for dw in [1,3,5]:
  p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5
  rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2;r=(rp+rn)/2
  scales[(dt,dw)]=np.maximum(1.4826*np.median(abs(r[15:-15]),axis=0),.002)
z=np.log1p(train);p=z-np.roll(z,1,axis=0);q=z-np.roll(z,-1,axis=0);r=((p-median_filter(p,size=(1,21)))+(q-median_filter(q,size=(1,21))))/2
scales['prior']=np.maximum(1.4826*np.median(abs(r[4:-4]),axis=0),.004)
med=median_filter(z,size=(9,1));scales['hampel']=np.maximum(1.4826*np.median(abs(z-med),axis=0),.002)
co=template_scores(train,scales);to=template_scores(test0,scales)
test,truth,events=inject(test0,20261126+sum(map(ord,h)),per=4);ts=template_scores(test,scales);cv=valid(train);tv=valid(test0)
rows=[];strata=[];real=[]
for model in co:
 th=float(np.quantile(co[model][cv],.999));pred=(ts[model]>th)&tv;orig=(to[model]>th)&tv
 evrec=lambda es:float(np.mean([pred[reg(e)].any() for e in es]))
 rows.append({'heat':h,'model':model,'threshold':th,'test_frames':len(test0),'events':len(events),'injected_points':int(truth.sum()),'point_recall':float(pred[truth].mean()),'event_recall':evrec(events),'original_candidate_count':int(orig.sum()),'candidate_rate':float(orig.sum()/tv.sum())})
 for dt,dw in shapes:
  for amp in [.05,.15,.4]:
   es=[e for e in events if e[2]==dt and e[3]==dw and e[5]==amp]
   strata.append({'heat':h,'model':model,'duration':dt,'width':dw,'amplitude':amp,'events':len(es),'event_recall':evrec(es)})
# peak proxy for multi local at 0.1% train threshold
th=float(np.quantile(co['multiscale_local'][cv],.999));flag=(to['multiscale_local']>th)&tv;cnt=hit=0
for i in range(15,len(test0)-15,5):
 p,_=find_peaks(test0[i],prominence=max(150,np.median(test0[i])*.08),width=(1,15),distance=5);p=p[(p>40)&(p<m-40)]
 for j in p:cnt+=1;hit+=int(flag[i,j])
real.append({'heat':h,'peak_proxy_n':cnt,'flagged_peak_proxy_n':hit,'flag_rate':hit/max(cnt,1)})
for name,arr in [('q1_anomaly_allheat.csv',rows),('q1_anomaly_strata.csv',strata),('q1_real_peak_proxy.csv',real)]:
 path=R/name;df=pd.DataFrame(arr)
 if path.exists():df=pd.concat([pd.read_csv(path),df],ignore_index=True)
 df.to_csv(path,index=False)
print(pd.DataFrame(rows).to_string(index=False));print('strata',pd.DataFrame(strata).groupby(['model','amplitude']).event_recall.mean().to_dict());print('peak proxy',real,flush=True)
