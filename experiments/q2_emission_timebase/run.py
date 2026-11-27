from pathlib import Path
import numpy as np, pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.special import expit
from scipy.optimize import least_squares
import argparse
parser=argparse.ArgumentParser(description="Q2.2 cumulative emission timebase sensitivity")
parser.add_argument('--input',type=Path,default=Path('data/interim/extracted/process.csv'))
parser.add_argument('--output',type=Path,default=Path('results/q2/emission_timebase/rerun'))
args=parser.parse_args(); root=args.output; root.mkdir(parents=True,exist_ok=True)
f=pd.read_csv(args.input,dtype={'heat':str}); D={}
assert sorted(f.attachment.unique())==list(range(1,11))
for h,frame in f.sort_values(['attachment','row']).groupby('heat',sort=False):
 assert frame.attachment.nunique()==1 and np.array_equal(frame.row.to_numpy(),np.arange(1,len(frame)+1))
 o=np.round(frame.O.to_numpy(dtype=float),2);y=np.round(frame[['CO2','CO']].to_numpy(dtype=float),2);D[h]=(o,y)
assert len(D)==10 and sum(len(o) for o,y in D.values())==13589 and list(D)[4]=='6C02057'
hs=list(D); gases=['CO2','CO']; grid=np.linspace(0,min(o[-1]-o[0] for o,y in D.values()),301)
def r2(y,p):return 1-np.sum((y-p)**2)/np.sum((y-y.mean())**2)
def fun(g,x,p):
 if g==0:
  a,b,v,w,u,s=p;return a*expit((x-v)/w)-b*expit((x-u)/s)
 c,v,w,u,s=p;return c*expit((x-v)/w)*expit((u-x)/s)
B=[([0,0,-500,1,-500,1],[120,120,4000,3000,4500,3000]),([0,-500,1,0,1],[150,4000,3000,5500,3000])]; P=[[27,14,124,160,950,400],[45,833,180,3294,120]]
def fit(train,g):
 x=np.concatenate([D[h][0][120::6] for h in train]);y=np.concatenate([D[h][1][120::6,g] for h in train]);w=np.concatenate([np.full(len(D[h][0][120::6]),1/np.sqrt(len(D[h][0][120::6]))) for h in train]);
 def res(p):
  z=fun(g,x,p);return np.r_[(z-y)*w,np.minimum(z,0)*w*3]
 return least_squares(res,P[g],bounds=B[g],max_nfev=100,ftol=1e-6).x
rows=[]; folds=[]; preds=[]
for g,gas in enumerate(gases):
 pp=fit(hs,g);lo={h:fit([j for j in hs if j!=h],g) for h in hs}
 for scen in ['H0','H1']:
  curves=[]
  for h in hs:
   o,y=D[h];dt=1 if scen=='H1' and h==hs[4] else .5
   x=np.maximum.accumulate(o-o[0]);e=cumulative_trapezoid(y[:,g],dx=dt,initial=0); # retain last point of each plateau
   k=np.r_[np.diff(x)>0,True];curves.append(np.interp(grid,x[k],e[k]))
  curves=np.array(curves)
  for mode in ['fit','LOHO']:
   pred=np.repeat(curves.mean(0)[None,:],10,axis=0) if mode=='fit' else (curves.sum(0)-curves)/9
   rr=[r2(y,p) for y,p in zip(curves,pred)]
   rows.append(dict(scenario=scen,gas=gas,model='mean_curve',mode=mode,macro_r2=np.mean(rr),pooled_r2=r2(curves.ravel(),pred.ravel())))
   for h,a,b,r in zip(hs,curves,pred,rr):
    folds.append(dict(scenario=scen,gas=gas,model='mean_curve',mode=mode,heat=h,r2=r));preds.extend(zip([scen]*len(grid),[gas]*len(grid),[mode]*len(grid),[h]*len(grid),grid,a,b))
  for mode in ['fit','LOHO']:
   obs=[];est=[];rr=[]
   for h in hs:
    o,y=D[h];dt=1 if scen=='H1' and h==hs[4] else .5
    # fixed row mask from prior H0 evaluation, not same physical time
    y=y[120:,g];v=np.maximum(fun(g,o[120:],pp if mode=='fit' else lo[h]),0)
    a=cumulative_trapezoid(y,dx=dt,initial=0);b=cumulative_trapezoid(v,dx=dt,initial=0);r=r2(a,b);rr.append(r);obs.extend(a);est.extend(b);folds.append(dict(scenario=scen,gas=gas,model='sigmoid_integral',mode=mode,heat=h,r2=r))
   rows.append(dict(scenario=scen,gas=gas,model='sigmoid_integral',mode=mode,macro_r2=np.mean(rr),pooled_r2=r2(np.array(obs),np.array(est))))
pd.DataFrame(rows).to_csv(root/'summary.csv',index=False);pd.DataFrame(folds).to_csv(root/'folds.csv',index=False);pd.DataFrame(preds,columns=['scenario','gas','mode','heat','oxygen_net','observed','prediction']).to_csv(root/'predictions.csv',index=False)
print('grid',grid[0],grid[-1]);print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(folds).query("model=='mean_curve' and mode=='LOHO'").pivot(index='heat',columns=['scenario','gas'],values='r2').to_string())
