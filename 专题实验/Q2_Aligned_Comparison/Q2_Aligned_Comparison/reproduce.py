from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.optimize import least_squares
from scipy.integrate import cumulative_trapezoid
R=Path(__file__).resolve().parent
D={p.stem[7:]:pd.read_csv(p) for p in sorted((R/'input').glob('curves_*.csv'))}
assert len(D)==10
G=['CO2','CO']; L=[0,12,24,36,48,60]
B={'CO2':([0,0,-500,1,-500,1],[120,120,4000,3000,4500,3000]),'CO':([0,-500,1,0,1],[150,4000,3000,5500,3000])}
P={'CO2':[27,14,124,160,950,400],'CO':[45,833,180,3294,120]}
def fun(g,x,p):
 if g=='CO2':
  a,b,v,w,u,s=p;return a*expit((x-v)/w)-b*expit((x-u)/s)
 c,v,w,u,s=p;return c*expit((x-v)/w)*expit((u-x)/s)
def get(h,lag,start=60):
 d=D[h];t=d.time_s.values; k=t>=start;return t[k],np.interp(t[k]-lag,t,d.oxygen_cumulative.values),{g:d[g+'_flow'].values[k] for g in G}
def fit(hs,lag,g,start=60):
 a=[get(h,lag,start) for h in hs];x=np.concatenate([v[1][::6] for v in a]);y=np.concatenate([v[2][g][::6] for v in a]);w=np.concatenate([np.full(len(v[1][::6]),1/np.sqrt(len(v[1][::6]))) for v in a]);
 def res(p):
  z=fun(g,x,p);return np.r_[(z-y)*w,np.minimum(z,0)*w*3]
 return least_squares(res,P[g],bounds=B[g],max_nfev=100,ftol=1e-6).x
def pred(g,x,p):return np.maximum(fun(g,x,p),0)
def met(y,p):
 y=np.asarray(y);p=np.asarray(p);e=y-p;return dict(r2=1-np.sum(e**2)/np.sum((y-y.mean())**2),mae=np.mean(abs(e)),rmse=np.sqrt(np.mean(e**2)))
rows=[];scan=[];pm=[];pars=[];predrows=[]
def plot(hs,name,lag,pp,start=60,held=None):
 fig,ax=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
 for j,g in enumerate(G):
  Y=[];Yp=[];E=[];Ep=[];ep=[];gx=np.linspace(max(get(h,lag,start)[1].min() for h in hs),min(get(h,lag,start)[1].max() for h in hs),300)
  for h in hs:
   t,x,y=get(h,lag,start);y=y[g];yp=pred(g,x,pp[g] if held is None else held[h][g]);e=cumulative_trapezoid(y,t,initial=0);e2=cumulative_trapezoid(yp,t,initial=0)
   line=ax[0,j].plot(x,y,lw=.9,alpha=.78,label=h)[0];ax[1,j].plot(x,e,lw=.9,alpha=.78)
   if held is not None:ax[0,j].plot(x,yp,'--',color=line.get_color(),lw=.85);ax[1,j].plot(x,e2,'--',color=line.get_color(),lw=.85)
   else:ep.append(np.interp(gx,x,e2));ax[1,j].plot(x,e2,'k--',alpha=.15,lw=.6)
   Y.extend(y);Yp.extend(yp);E.extend(e);Ep.extend(e2)
  m=met(Y,Yp);me=met(E,Ep)
  if held is None:
   xx=np.linspace(min(get(h,lag,start)[1].min() for h in hs),max(get(h,lag,start)[1].max() for h in hs),500);ax[0,j].plot(xx,pred(g,xx,pp[g]),'k--',lw=2,label='Pooled sigmoid');ax[1,j].plot(gx,np.mean(ep,axis=0),'k--',lw=2,label='Mean integrated model');ax[1,j].legend(fontsize=8)
  ax[0,j].legend(fontsize=7,ncol=2);ax[0,j].set_title(f'{g} flow | pooled R²={m["r2"]:.3f}');ax[1,j].set_title(f'{g} integrated volume | pooled R²={me["r2"]:.3f}');ax[0,j].set_ylabel('Volumetric flow (m³/s)');ax[1,j].set_ylabel(f'Volume since t={start} s (m³)')
  for a in ax[:,j]:a.set_xlabel(f'Cumulative O₂ at t−{lag} s (m³)');a.grid(alpha=.2)
  pm.extend([dict(figure=name,gas=g,quantity='flow',**m),dict(figure=name,gas=g,quantity='integrated',**me)])
 fig.suptitle(f'{name} | '+('held-out heat: solid observed / dashed predicted' if held else 'descriptive in-sample fit'),fontsize=13);fig.savefig(R/f'{name}.png',dpi=190);plt.close(fig)
for name,hs in [('all10',list(D)),('att2to10',[h for h in D if h!='6B03205'])]:
 for lag in L:
  pp={g:fit(hs,lag,g) for g in G}
  for g in G:
   m=[met(get(h,lag)[2][g],pred(g,get(h,lag)[1],pp[g])) for h in hs];scan.append(dict(group=name,lag_s=lag,gas=g,mean_heat_r2=np.mean([a['r2'] for a in m]),mean_heat_mae=np.mean([a['mae'] for a in m])))
  if lag in [0,48]:
   plot(hs,f'{name}_lag{lag:02d}_four_panel',lag,pp)
   for g in G:
    pars.append(dict(group=name,lag_s=lag,gas=g,parameters=';'.join(map(str,pp[g]))))
    for h in hs:
     t,x,y=get(h,lag);y=y[g];pi=fit([h],lag,g)
     for mode,pr in [('pooled_in_sample',pred(g,x,pp[g])),('per_heat_in_sample',pred(g,x,pi))]:rows.append(dict(group=name,lag_s=lag,mode=mode,heat=h,gas=g,**met(y,pr)))
 if name=='att2to10':
  plot(hs,'att2to10_lag48_full_aligned_domain',48,{g:fit(hs,48,g,48) for g in G},start=48)
 held={}
 for h in hs:
  tr=[a for a in hs if a!=h];cand=[]
  for lag in L:
   pp={g:fit(tr,lag,g) for g in G};score=sum(np.mean([1-met(get(a,lag)[2][g],pred(g,get(a,lag)[1],pp[g]))['r2'] for a in tr]) for g in G);cand.append((score,lag,pp))
  chosen=min(cand,key=lambda a:a[0]);held[h]=next(a[2] for a in cand if a[1]==48)
  for mode,c in [('LOHO_fixed0',cand[0]),('LOHO_fixed48',next(a for a in cand if a[1]==48)),('LOHO_train_selected',chosen)]:
   _,lag,pp=c;t,x,y=get(h,lag)
   for g in G:
    yp=pred(g,x,pp[g]);rows.append(dict(group=name,lag_s=lag,mode=mode,heat=h,gas=g,**met(y[g],yp)));predrows.extend(zip([name]*len(t),[mode]*len(t),[h]*len(t),[g]*len(t),[lag]*len(t),t,x,y[g],yp))
 plot(hs,f'{name}_lag48_LOHO_four_panel',48,None,held=held)
 print(name,'done',flush=True)
A=pd.DataFrame(rows);A.to_csv(R/'metrics_by_heat.csv',index=False);pd.DataFrame(scan).to_csv(R/'lag_scan_descriptive.csv',index=False);pd.DataFrame(pm).to_csv(R/'pooled_figure_metrics.csv',index=False);pd.DataFrame(pars).to_csv(R/'pooled_parameters.csv',index=False)
A[~A['mode'].str.startswith('LOHO')].groupby(['group','mode','gas','lag_s'])[['r2','mae','rmse']].mean().to_csv(R/'in_sample_summary.csv');A[A['mode'].str.startswith('LOHO')].groupby(['group','mode','gas'])[['r2','mae','rmse']].mean().to_csv(R/'LOHO_summary.csv');pd.DataFrame(predrows,columns=['group','mode','heat','gas','lag_s','time_s','oxygen_aligned','observed','predicted']).to_csv(R/'LOHO_predictions.csv',index=False)
print(pd.read_csv(R/'in_sample_summary.csv').to_string(index=False));print(pd.read_csv(R/'LOHO_summary.csv').to_string(index=False))
