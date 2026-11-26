import os,glob
import numpy as np,pandas as pd
from scipy.interpolate import BSpline,PchipInterpolator,UnivariateSpline
from scipy.optimize import lsq_linear
from scipy.stats import spearmanr
R='/mnt/data/APMCM_Full_Experiments';files=sorted(glob.glob('/mnt/data/heats_*.npz'))
D={}
for f in files:
 z=np.load(f);h=os.path.basename(f)[6:-4];D[h]={k:z[k] for k in z.files}
def basis_int(t,spacing):
 knots=np.r_[[t[0]]*4,np.arange(t[0]+spacing,t[-1],spacing),[t[-1]]*4];p=len(knots)-4
 bs=BSpline(knots,np.eye(p),3,extrapolate=True);B=np.asarray(bs(t));anti=bs.antiderivative();A=np.asarray(anti(t)-anti(t[0]));return B,A,p
def fitflow(O,t,spacing,lam,mask):
 B,A,p=basis_int(t,spacing);D2=np.diff(np.eye(p),n=2,axis=0);w=np.ones(mask.sum());coef=None;ok=False
 for _ in range(2):
  M=np.vstack([A[mask]*np.sqrt(w)[:,None],np.sqrt(lam)*D2]);y=np.r_[(O[mask]-O[0])*np.sqrt(w),np.zeros(len(D2))]
  fit=lsq_linear(M,y,bounds=(0,np.inf),max_iter=80,tol=1e-7);coef=fit.x;ok=fit.success
  r=O[mask]-O[0]-A[mask]@coef;s=max(1.4826*np.median(abs(r-np.median(r))),.1);w=np.minimum(1,1.345*s/np.maximum(abs(r),1e-9))
 return B@coef,O[0]+A@coef,ok

def blocks(n):
 out=[]
 for frac in [.12,.27,.42,.57,.72,.87]:
  a=max(20,int(n*frac));b=min(n-20,a+16);out.append(np.arange(a,b))
 return out
flowcv=[];gascv=[];supports=[];counter=[];phases=[];assoc=[];curvefiles=[]
for h,z in D.items():
 n=len(z['seq']);t=np.arange(n)*.5;O=z['proc'][:,0].astype(float);G=z['proc'][:,1:3].astype(float);bs=blocks(n);outer=np.r_[bs[1],bs[4]];outermask=np.ones(n,bool);outermask[outer]=False;inner=[bs[k] for k in [0,2,3,5]]
 # inner tune lambda and knot spacing, excluding outer test blocks from all inner fits
 candidates=[]
 for spacing in [5,10,20]:
  for lam in [.1,1,10,100]:
   es=[]
   for ix in inner:
    tr=outermask.copy();tr[ix]=False
    _,oh,_=fitflow(O,t,spacing,lam,tr);es.extend(oh[ix]-O[ix])
   candidates.append((np.mean(np.abs(np.asarray(es))),spacing,lam))
 _,spacing,lam=min(candidates,key=lambda a:a[0])
 q,oh,ok=fitflow(O,t,spacing,lam,outermask)
 # honest outer test for cumulative oxygen
 lin=np.interp(t[outer],t[outermask],O[outermask]);pch=PchipInterpolator(t[outermask],O[outermask])(t[outer])
 qcv,ohcv,okcv=fitflow(O,t,spacing,lam,outermask)
 for model,pred in [('linear',lin),('PCHIP',pch),('nonnegative_integral_spline',ohcv[outer])]:
  flowcv.append({'heat':h,'model':model,'outer_MAE':float(np.mean(abs(pred-O[outer]))),'outer_RMSE':float(np.sqrt(np.mean((pred-O[outer])**2))),'test_points':len(outer),'spacing_s':spacing if model.endswith('spline') else np.nan,'penalty':lam if model.endswith('spline') else np.nan})
 # full data q derivatives and support
 q_diff=np.r_[np.diff(O)/.5,np.nan];q_p=np.maximum(PchipInterpolator(t,O).derivative()(t),0)
 def slope(win,causal=False):
  k=int(win/.5);out=np.full(n,np.nan)
  for i in range(n):
   a=max(0,i-k+1) if causal else max(0,i-k//2);b=i+1 if causal else min(n,i+k//2+1)
   if b-a<3:continue
   x=t[a:b];y=O[a:b];xc=x-x.mean();out[i]=max(0,xc@(y-y.mean())/(xc@xc))
  return out
 q2=slope(2);q5=slope(5);q10=slope(10);qtrail=slope(10,True)
 for name,qq in [('difference',q_diff),('PCHIP_derivative',q_p),('integral_spline',q),('centered_slope_2s',q2),('centered_slope_5s',q5),('centered_slope_10s',q10),('trailing_slope_10s',qtrail)]:
  v=np.isfinite(qq);supports.append({'heat':h,'model':name,'p10':np.nanquantile(qq[v],.1),'p50':np.nanquantile(qq[v],.5),'p90':np.nanquantile(qq[v],.9),'sd':np.nanstd(qq[v]),'negative_count':int(np.sum(qq[v]<-1e-8))})
 dO=np.diff(O);upd=np.flatnonzero(dO>0)+1;gaps=np.diff(np.r_[0,upd])*.5;pos=dO[dO>0];med=np.median(pos) if len(pos) else np.nan
 counter.append({'heat':h,'flat_fraction':np.mean(dO==0),'positive_updates':len(upd),'median_increment':med,'median_update_gap_s':np.median(gaps) if len(gaps) else np.nan,'update_gap_1s_share':np.mean(np.isclose(gaps,1)) if len(gaps) else np.nan,'near_integer_median_increment_share':np.mean(abs(pos/med-np.round(pos/med))<.04) if len(pos) and med else np.nan})
 # Gas nested block CV, using same outer blocks for all methods
 for j,gn in enumerate(['CO2','CO']):
  y=G[:,j];opts=[]
  for s in [0,n*.01,n*.1,n,10*n]:
   es=[]
   for ix in inner:
    tr=outermask.copy();tr[ix]=False;spl=UnivariateSpline(t[tr],y[tr],s=s,k=3);es.extend(spl(t[ix])-y[ix])
   opts.append((np.mean(np.abs(np.asarray(es))),s))
  _,best_s=min(opts,key=lambda a:a[0])
  # fit on non-outer, evaluate on outer
  p_lin=np.interp(t[outer],t[outermask],y[outermask]);p_pch=PchipInterpolator(t[outermask],y[outermask])(t[outer]);p_spl=np.maximum(0,UnivariateSpline(t[outermask],y[outermask],s=best_s,k=3)(t[outer]))
  for model,pred in [('linear',p_lin),('PCHIP',p_pch),('smoothing_spline',p_spl)]:
   gascv.append({'heat':h,'gas':gn,'model':model,'outer_MAE':float(np.mean(abs(pred-y[outer]))),'outer_RMSE':float(np.sqrt(np.mean((pred-y[outer])**2))),'n_test':len(outer),'smoothing_s':best_s if model=='smoothing_spline' else np.nan})
 # phase thirds, cumulative gas and associations
 c2=np.r_[0,np.cumsum((G[:-1,0]+G[1:,0])*.25)];c1=np.r_[0,np.cumsum((G[:-1,1]+G[1:,1])*.25)];tot=G.sum(1);share=np.divide(G[:,0],tot,out=np.full(n,np.nan),where=tot>.1);eq=G[:,0]+G[:,1]/2
 for k in range(3):
  a=k*n//3;b=(k+1)*n//3;sl=slice(a,b);phases.append({'heat':h,'time_third':k+1,'t_start_s':t[a],'t_end_s':t[b-1],'q10_median':np.nanmedian(q10[sl]),'q10_sd':np.nanstd(q10[sl]),'CO2_median':np.median(G[sl,0]),'CO_median':np.median(G[sl,1]),'CO2_share_median':np.nanmedian(share[sl]),'CO2_volume_m3':c2[b-1]-c2[a],'CO_volume_m3':c1[b-1]-c1[a]})
 for name,y in [('CO2',G[:,0]),('CO',G[:,1]),('total',tot),('CO2_share',share),('oxygen_equivalent',eq)]:
  v=np.isfinite(y);assoc.append({'heat':h,'variable':name,'spearman_O_level':spearmanr(O[v],y[v]).statistic,'spearman_dO_dy':spearmanr(np.diff(O[v]),np.diff(y[v])).statistic,'n':int(v.sum())})
 # save continuous curves from selected nonnegative flow model and raw gas signals
 pd.DataFrame({'time_s':t,'oxygen_cumulative':O,'oxygen_flow_spline':q,'oxygen_flow_centered10':q10,'oxygen_flow_trailing10':qtrail,'CO2_flow':G[:,0],'CO_flow':G[:,1],'CO2_cumulative':c2,'CO_cumulative':c1,'CO2_plus_CO':tot,'CO2_share':share}).to_csv(R+f'/curves_{h}.csv',index=False)
 print('nested CV complete',h,'spacing',spacing,'lambda',lam,flush=True)
for name,rows in [('q2_flow_nested_cv.csv',flowcv),('q2_flow_support_nested.csv',supports),('q2_counter_updates_nested.csv',counter),('q2_gas_nested_cv.csv',gascv),('q2_phase_nested.csv',phases),('q2_associations_nested.csv',assoc)]:pd.DataFrame(rows).to_csv(R+'/'+name,index=False)
print('\nFLOW NESTED\n',pd.DataFrame(flowcv).groupby('model')[['outer_MAE','outer_RMSE']].mean().round(4).to_string())
print('\nGAS NESTED\n',pd.DataFrame(gascv).groupby(['gas','model'])[['outer_MAE','outer_RMSE']].mean().round(4).to_string())
