from pathlib import Path
import os,glob,json
import numpy as np,pandas as pd
from scipy.ndimage import uniform_filter1d as uf, maximum_filter,median_filter,shift as ndi_shift
from scipy.signal import find_peaks
R=Path('/mnt/data/APMCM_Full_Experiments')
shapes=[(1,1),(1,3),(1,5),(3,1),(3,3),(5,5)]
def region(e):
 i,j,dt,dw,*_=e;return slice(i-dt//2,i+dt//2+1),slice(j-dw//2,j+dw//2+1)
def inject(a,seed,per=8):
 rg=np.random.default_rng(seed);b=a.copy();y=np.zeros(a.shape,bool);ev=[];m=a.shape[1]
 for dt,dw in shapes:
  for s in [1,-1]:
   for amp in [.05,.15,.4]:
    for _ in range(per):
     for trial in range(10000):
      i=int(rg.integers(20,len(a)-20));j=int(rg.integers(45,m-45));e=(i,j,dt,dw,s,amp);ss=region(e)
      if y[i-16:i+17,j-35:j+36].any() or a[ss].min()<300 or a[ss].max()>45000:continue
      break
     else:continue
     b[ss]=a[ss]*(1+s*amp);y[ss]=True;ev.append(e)
 return b,y,ev
def valid(a):
 v=np.ones(a.shape,bool);v[:15]=v[-15:]=False;v[:,:40]=v[:,-40:]=False;return v
def templates(a):
 z=np.log1p(a)
 for dt in [1,3,5]:
  za=uf(z,dt,axis=0,mode='nearest');p=za-np.roll(za,dt,axis=0);q=za-np.roll(za,-dt,axis=0)
  for dw in [1,3,5]:
   p1=uf(p,dw,axis=1);q1=uf(q,dw,axis=1);pb=uf(p,9,axis=1);qb=uf(q,9,axis=1);off=dw//2+5
   rp=p1-(np.roll(pb,off,axis=1)+np.roll(pb,-off,axis=1))/2;rn=q1-(np.roll(qb,off,axis=1)+np.roll(qb,-off,axis=1))/2
   r=(rp+rn)/2;ag=np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0);yield (dt,dw),r,ag
def prior(a):
 z=np.log1p(a);p=z-np.roll(z,1,axis=0);q=z-np.roll(z,-1,axis=0);rp=p-median_filter(p,size=(1,21));rn=q-median_filter(q,size=(1,21));return (rp+rn)/2,np.where(rp*rn>0,np.minimum(abs(rp),abs(rn)),0)
def make_scores(train,test):
 scales={sh:np.maximum(1.4826*np.median(abs(r[15:-15]),axis=0),.002) for sh,r,ag in templates(train)}
 r,_=prior(train);scales['prior']=np.maximum(1.4826*np.median(abs(r[4:-4]),axis=0),.004)
 def scores(a):
  out={k:np.zeros(a.shape,np.float32) for k in ['single_local','multiscale','multiscale_local','prior_soft']}
  for sh,r,ag in templates(a):
   dt,dw=sh;v=.65*abs(r)+.35*ag;base=scales[sh];s=v/base;w=uf(abs(r),41,axis=1);off=dw//2+21;local=(np.roll(w,off,axis=1)+np.roll(w,-off,axis=1))/2;sl=v/np.maximum(base,1.4826*local)
   if sh==(1,1):out['single_local']=sl.copy()
   np.maximum(out['multiscale'],maximum_filter(s,size=sh),out=out['multiscale']);np.maximum(out['multiscale_local'],maximum_filter(sl,size=sh),out=out['multiscale_local'])
  r,ag=prior(a);out['prior_soft']=(.65*abs(r)+.35*ag)/np.maximum(np.maximum(1.4826*median_filter(abs(r),size=(1,21)),scales['prior']),.004)
  return out
 return scores(train),scores(test)
def event_recall(mask,events):return float(np.mean([mask[region(e)].any() for e in events])) if events else np.nan
allrows=[]; strata=[]; realrows=[]
for f in sorted(glob.glob('/mnt/data/heats_*.npz')):
 z=np.load(f);h=os.path.basename(f)[6:-4];X=z['X'].astype(np.float32);n=len(X);cut=int(.55*n);gap=40;train=X[:cut];test0=X[min(cut+gap,n-30):]
 if len(test0)<200: test0=X[cut:]
 test,truth,events=inject(test0,20261126+len(allrows),per=8)
 co,ct=make_scores(train,test)
 cv=valid(train);tv=valid(test0)
 # Add Hampel temporal median score
 med=median_filter(np.log1p(train),size=(9,1));sc=np.maximum(1.4826*np.median(abs(np.log1p(train)-med),axis=0),.002)
 def hamp(a):return abs(np.log1p(a)-median_filter(np.log1p(a),size=(9,1)))/sc
 co['hampel']=hamp(train);ct['hampel']=hamp(test)
 for model in co:
  budget=.001;th=float(np.quantile(co[model][cv],1-budget));pred=(ct[model]>th)&tv
  allrows.append({'heat':h,'model':model,'threshold':th,'injected_points':int(truth.sum()),'point_recall':float(pred[truth].mean()),'event_recall':event_recall(pred,events),'original_flags':int(((co[model]>th)&cv).sum()),'test_original_flags':int(((make_scores(train,test0)[1][model]>th)&tv).sum()),'test_candidate_rate':float(((make_scores(train,test0)[1][model]>th)&tv).sum()/tv.sum())})
  for dt,dw in shapes:
   for amp in [.05,.15,.4]:
    es=[e for e in events if e[2]==dt and e[3]==dw and e[5]==amp]
    if es:
     strata.append({'heat':h,'model':model,'duration':dt,'width':dw,'amplitude':amp,'events':len(es),'event_recall':event_recall(pred,es),'point_recall':float(pred[truth & np.isin(np.indices(truth.shape)[0],[])].mean()) if False else np.nan})
 # proxy of repeated real peaks: count candidates flagged at a detected peak in 1-frame neighborhoods
 scoretest=make_scores(train,test0)[1]
 mname='multiscale_local';th=float(np.quantile(co[mname][cv],.999));
 flagged=(scoretest[mname]>th)&tv;cnt=0;hit=0
 for i in range(15,len(test0)-15,5):
  p,_=find_peaks(test0[i],prominence=max(150,np.median(test0[i])*.08),width=(1,15),distance=5)
  p=p[(p>40)&(p<test0.shape[1]-40)]
  for j in p:
   cnt+=1;hit+=int(flagged[i,j])
 realrows.append({'heat':h,'persistent_peak_proxy_n':cnt,'flagged_peak_proxy_n':hit,'proxy_flag_rate':hit/max(cnt,1)})
 print('done',h,'rows',n,'events',len(events),flush=True)
pd.DataFrame(allrows).to_csv(R/'q1_anomaly_allheat.csv',index=False)
pd.DataFrame(strata).to_csv(R/'q1_anomaly_strata.csv',index=False)
pd.DataFrame(realrows).to_csv(R/'q1_real_peak_proxy.csv',index=False)
print('\nSUMMARY')
print(pd.DataFrame(allrows).groupby('model').agg(heat_mean_point_recall=('point_recall','mean'),heat_mean_event_recall=('event_recall','mean'),heat_sd_event_recall=('event_recall','std'),mean_test_candidate_rate=('test_candidate_rate','mean')).round(4).to_string())
print('\nPER HEAT')
print(pd.DataFrame(allrows).pivot(index='heat',columns='model',values='event_recall').round(3).to_string())
print('\nSTRATA')
print(pd.DataFrame(strata).groupby(['model','amplitude']).event_recall.mean().unstack().round(3).to_string())
