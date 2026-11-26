"""Shared numerical helpers for the APMCM time-base audit."""
import math
import numpy as np
from scipy.signal import savgol_filter, find_peaks, welch, correlate
from scipy.interpolate import BSpline
from scipy.optimize import lsq_linear
from scipy.ndimage import uniform_filter1d


def acf(x,maxlag=200):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 if len(x)<5: return np.full(maxlag+1,np.nan)
 x=x-np.mean(x); den=x@x
 if den<=0: return np.full(maxlag+1,np.nan)
 out=np.full(maxlag+1,np.nan); out[0]=1.
 for k in range(1,min(maxlag+1,len(x))): out[k]=x[:-k]@x[k:]/den
 return out

def decay_lag(a,threshold=1/math.e):
 a=np.asarray(a); ix=np.where(np.isfinite(a[1:])&(a[1:]<=threshold))[0]
 return int(ix[0]+1) if len(ix) else np.nan

def psd_summary(x,dt):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 if len(x)<16 or np.std(x)<1e-12: return {'psd_peak_hz':np.nan,'psd_median_hz':np.nan,'psd_centroid_hz':np.nan}
 f,p=welch(x-np.mean(x),fs=1/dt,nperseg=min(256,len(x)))
 if len(f)<2 or p[1:].sum()<=0: return {'psd_peak_hz':np.nan,'psd_median_hz':np.nan,'psd_centroid_hz':np.nan}
 ff=f[1:]; pp=p[1:]; c=np.cumsum(pp); med=ff[np.searchsorted(c,c[-1]/2)]
 return {'psd_peak_hz':float(ff[np.argmax(pp)]),'psd_median_hz':float(med),'psd_centroid_hz':float((ff@pp)/pp.sum())}

def run_lengths(mask):
 m=np.asarray(mask,bool)
 if not m.any(): return []
 d=np.diff(np.r_[False,m,False].astype(int))
 return (np.where(d==-1)[0]-np.where(d==1)[0]).tolist()

def event_duration(x,dt,quantile=.90):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 if len(x)<5: return {'event_count':0,'event_duration_median_s':np.nan,'event_duration_p90_s':np.nan}
 ll=run_lengths(x>=np.quantile(x,quantile))
 if not ll: return {'event_count':0,'event_duration_median_s':0.,'event_duration_p90_s':0.}
 a=np.asarray(ll)*dt
 return {'event_count':len(a),'event_duration_median_s':float(np.median(a)),'event_duration_p90_s':float(np.quantile(a,.9))}

def spectral_features(X,wl):
 X=np.asarray(X,float); wl=np.asarray(wl,float); sums=X.sum(axis=1)
 integ=np.trapezoid(X,wl,axis=1); mean=X.mean(axis=1); p95=np.quantile(X,.95,axis=1); mx=X.max(axis=1)
 centroid=np.divide(X@wl,sums,out=np.full(len(X),np.nan),where=sums>0); sat=np.mean(X==65535,axis=1)
 a=X[:-1]; b=X[1:]; den=np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)
 cosine=np.divide(np.sum(a*b,axis=1),den,out=np.full(len(a),np.nan),where=den>0)
 ac=a-a.mean(axis=1,keepdims=True); bc=b-b.mean(axis=1,keepdims=True); den2=np.linalg.norm(ac,axis=1)*np.linalg.norm(bc,axis=1)
 corr=np.divide(np.sum(ac*bc,axis=1),den2,out=np.full(len(a),np.nan),where=den2>0)
 med=np.median(X,axis=0); sm=savgol_filter(med,15,2,mode='interp'); prom=max(float(np.percentile(sm,95)-np.percentile(sm,5))*.025,1.)
 peaks,props=find_peaks(sm,prominence=prom,distance=12)
 if len(peaks): peaks=peaks[np.argsort(props['prominences'])[::-1][:5]]
 feat={'spectral_integral':integ,'spectral_mean':mean,'spectral_p95':p95,'spectral_max':mx,'spectral_centroid':centroid,
       'saturation_fraction':sat,'adjacent_cosine':np.r_[cosine,np.nan],'adjacent_pearson':np.r_[corr,np.nan]}
 for j,p in enumerate(peaks,1):
  lo=max(0,int(p)-18); hi=min(X.shape[1],int(p)+19); xx=wl[lo:hi]; block=X[:,lo:hi]
  loc=np.argmax(block,axis=1); h=block[np.arange(len(X)),loc]; pos=xx[loc]; base=np.minimum(block[:,0],block[:,-1])
  width=np.zeros(len(X)); area=np.zeros(len(X))
  for i in range(len(X)):
   half=base[i]+.5*(h[i]-base[i]); width[i]=np.sum(block[i]>=half)*np.median(np.diff(xx))
   area[i]=np.trapezoid(np.maximum(block[i]-base[i],0),xx)
  feat[f'peak{j}_height']=h; feat[f'peak{j}_position_nm']=pos; feat[f'peak{j}_width_nm']=width; feat[f'peak{j}_area']=area
 return feat

def flow_slope(O,dt,win_s=10.,causal=False):
 O=np.asarray(O,float); n=len(O); k=max(3,int(round(win_s/dt))); left=k//2; right=k//2 if k%2==0 else k//2+1
 q=np.full(n,np.nan)
 for i in range(n):
  a=max(0,i-k+1) if causal else max(0,i-left); b=i+1 if causal else min(n,i+right)
  if b-a<3: continue
  x=np.arange(a,b,dtype=float)*dt; y=O[a:b]; xc=x-x.mean(); den=xc@xc
  if den>0: q[i]=max(0.,float(xc@(y-y.mean())/den))
 return q

def integral_basis(t,spacing):
 t=np.asarray(t,float); inner=np.arange(t[0]+spacing,t[-1],spacing); knots=np.r_[[t[0]]*4,inner,[t[-1]]*4]
 p=len(knots)-4; bs=BSpline(knots,np.eye(p),3,extrapolate=True); B=np.asarray(bs(t)); anti=bs.antiderivative(); A=np.asarray(anti(t)-anti(t[0]))
 return B,A,p

def fit_integral_flow(O,t,spacing,lam,trainmask):
 O=np.asarray(O,float); B,A,p=integral_basis(t,spacing); D2=np.diff(np.eye(p),n=2,axis=0); w=np.ones(int(trainmask.sum())); coef=None
 for _ in range(3):
  M=np.vstack([A[trainmask]*np.sqrt(w)[:,None],np.sqrt(lam)*D2]); y=np.r_[(O[trainmask]-O[0])*np.sqrt(w),np.zeros(len(D2))]
  fit=lsq_linear(M,y,bounds=(0,np.inf),max_iter=100,tol=1e-8); coef=fit.x; r=O[trainmask]-O[0]-A[trainmask]@coef
  s=max(1.4826*np.median(np.abs(r-np.median(r))),.05); w=np.minimum(1.,1.345*s/np.maximum(np.abs(r),1e-9))
 return B@coef,O[0]+A@coef

def block_sets(n):
 blocks=[]
 for frac in [.12,.27,.42,.57,.72,.87]:
  a=max(20,int(n*frac)); b=min(n-20,a+16)
  if b>a: blocks.append(np.arange(a,b))
 return blocks

def fit_flow_cv(O,dt):
 n=len(O); t=np.arange(n)*dt; blocks=block_sets(n); outer=np.r_[blocks[1],blocks[4]]; inner=[blocks[k] for k in [0,2,3,5]]
 mask=np.ones(n,bool); mask[outer]=False; candidates=[]
 for spacing in [5.,10.,20.]:
  for lam in [.1,1.,10.,100.]:
   e=[]
   for ix in inner:
    tr=mask.copy(); tr[ix]=False; _,pred=fit_integral_flow(O,t,spacing,lam,tr); e.extend(pred[ix]-O[ix])
   candidates.append((float(np.mean(np.abs(e))),spacing,lam))
 best=min(candidates); _,spacing,lam=best; q,oh=fit_integral_flow(O,t,spacing,lam,mask)
 lin=np.interp(t[outer],t[mask],O[mask]); pch=__import__('scipy').interpolate.PchipInterpolator(t[mask],O[mask])(t[outer])
 full=np.ones(n,bool); qfull,ofull=fit_integral_flow(O,t,spacing,lam,full)
 return qfull,ofull,{'dt_s':dt,'n_test':len(outer),'spacing_s':spacing,'lambda':lam,'inner_MAE':best[0],
  'integral_spline_outer_MAE':float(np.mean(np.abs(oh[outer]-O[outer]))),'integral_spline_outer_RMSE':float(np.sqrt(np.mean((oh[outer]-O[outer])**2))),
  'linear_outer_MAE':float(np.mean(np.abs(lin-O[outer]))),'pchip_outer_MAE':float(np.mean(np.abs(pch-O[outer])))}

def corr_lags(x,y,maxlag):
 x=np.asarray(x,float); y=np.asarray(y,float); n=min(len(x),len(y)); x=x[:n]; y=y[:n]
 mx=np.isfinite(x).astype(float); my=np.isfinite(y).astype(float); xx=np.nan_to_num(x); yy=np.nan_to_num(y)
 def c(a,b): return correlate(a,b,mode='full',method='fft')[n-1-maxlag:n][::-1]
 count=c(mx,my); sx=c(xx*mx,my); sy=c(mx,yy*my); sxx=c(xx*xx*mx,my); syy=c(mx,yy*yy*my); sxy=c(xx*mx,yy*my)
 with np.errstate(divide='ignore',invalid='ignore'):
  cov=sxy-sx*sy/count; vx=sxx-sx*sx/count; vy=syy-sy*sy/count
  out=cov/np.sqrt(np.maximum(vx,0)*np.maximum(vy,0))
 out[(count<50)|(vx<=1e-12)|(vy<=1e-12)]=np.nan
 return out

def lag_estimate(q,gas,dt,seed=20261126):
 win=max(3,int(round(5/dt))); qs=uniform_filter1d(np.nan_to_num(q,nan=np.nanmedian(q)),size=win,mode='nearest'); gs=uniform_filter1d(np.asarray(gas,float),size=win,mode='nearest')
 dq=np.diff(qs,prepend=qs[0]); dy=np.diff(gs,prepend=gs[0]); maxlag=max(2,int(round(40/dt))); cc=corr_lags(dq,dy,maxlag)
 if not np.isfinite(cc).any(): return {'lag_s':np.nan,'lag_rows':np.nan,'corr':np.nan,'boundary':True,'shift_p':np.nan,'lag0_corr':np.nan,'max_abs_corr':np.nan}
 k=int(np.nanargmax(np.abs(cc))); obs=float(np.nanmax(np.abs(cc))); rng=np.random.default_rng(seed); null=[]
 lo=maxlag+2; hi=len(dy)-maxlag-2
 if hi>lo:
  for s in rng.integers(lo,hi,size=99):
   c0=corr_lags(dq,np.roll(dy,int(s)),maxlag)
   if np.isfinite(c0).any(): null.append(np.nanmax(np.abs(c0)))
 p=(1+sum(v>=obs for v in null))/(1+len(null)) if null else np.nan
 return {'lag_s':k*dt,'lag_rows':k,'corr':float(cc[k]),'boundary':bool(k in [0,maxlag]),'shift_p':float(p),
  'lag0_corr':float(cc[0]),'max_abs_corr':obs}
