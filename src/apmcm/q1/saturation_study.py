import os,glob
import numpy as np,pandas as pd
from scipy.ndimage import label
from scipy.signal import find_peaks
from scipy.optimize import minimize
from scipy.special import log_ndtr
from apmcm.common.paths import PROC, RESULTS
from apmcm.q1.reconstruct import SATURATION_CEILING
R=RESULTS/'q1/rerun/saturation_study/tables'; R.mkdir(parents=True,exist_ok=True)
files=sorted(str(p) for p in PROC.glob('raw_*.npz'))
def run():
    summary=[];regions=[];simrows=[]
    for f in files:
     z=np.load(f);h=os.path.basename(f)[4:-4];X=z['x'].astype(np.float32);wl=z['wl'];n,m=X.shape;sat=X>=SATURATION_CEILING
     if sat.any():
      # per-pixel temporal availability
      prev=np.zeros_like(sat);nxt=np.zeros_like(sat);prev[1:]=~sat[:-1];nxt[:-1]=~sat[1:]
      both=sat&prev&nxt;one=sat&(prev|nxt)&~both;none=sat&~(prev|nxt)
      lab,nc=label(sat,structure=np.ones((3,3),int));sizes=np.bincount(lab.ravel())[1:];
      for k in range(1,nc+1):
       ii,jj=np.where(lab==k);regions.append({'heat':h,'component':k,'size':len(ii),'t_start_s':ii.min()*.5,'t_end_s':ii.max()*.5,'lambda_min':wl[jj.min()],'lambda_max':wl[jj.max()],'time_span_s':(ii.max()-ii.min())*.5,'spectral_span_nm':wl[jj.max()]-wl[jj.min()]})
      summary.append({'heat':h,'sat_points':int(sat.sum()),'sat_rows':int(sat.any(1).sum()),'sat_columns':int(sat.any(0).sum()),'both_adjacent_frames_unsat':int(both.sum()),'one_adjacent_frame_unsat':int(one.sum()),'neither_adjacent_frame_unsat':int(none.sum()),'both_share':float(both.sum()/sat.sum()),'one_share':float(one.sum()/sat.sum()),'neither_share':float(none.sum()/sat.sum()),'components':nc,'largest_component':int(sizes.max())})
     else:
      summary.append({'heat':h,'sat_points':0,'sat_rows':0,'sat_columns':0,'both_adjacent_frames_unsat':0,'one_adjacent_frame_unsat':0,'neither_adjacent_frame_unsat':0,'both_share':np.nan,'one_share':np.nan,'neither_share':np.nan,'components':0,'largest_component':0})
     # choose unsaturated peak candidates for artificial censoring
     rng=np.random.default_rng(20261126+sum(map(ord,h)));cand=[]
     for i in range(20,n-20,10):
      if sat[i].any() or sat[i-1].any() or sat[i+1].any():continue
      p,props=find_peaks(X[i],prominence=max(250,np.median(X[i])*.12),width=(1,15),distance=5)
      for j in p:
       if j<30 or j>=m-30:continue
       if np.any(X[i,j-12:j+13]>=65535) or np.any(X[i-1,j-12:j+13]>=65535) or np.any(X[i+1,j-12:j+13]>=65535):continue
       cand.append((i,int(j)))
     rng.shuffle(cand);cand=cand[:20]
     for i,j in cand:
      xx=wl[j-12:j+13]-wl[j];a=X[i,j-12:j+13].astype(float);base=np.linspace(a[0],a[-1],len(a));hgt=max(a.max()-base[12],1)
      shoulders=np.r_[0:6,19:25];sig=max(5,1.4826*np.median(abs(np.diff(a[shoulders])-np.median(np.diff(a[shoulders]))))/np.sqrt(2))
      for frac in [.5,.75]:
       cap=base[12]+frac*hgt;obs=np.minimum(a,cap);mask=a>=cap
       if mask.sum()==0:continue
       # linear interpolation across censored interval
       predlin=np.interp(np.arange(len(a)),np.flatnonzero(~mask),obs[~mask])
       # local shoulder-corrected temporal template
       sh=np.r_[np.arange(3,10),np.arange(15,22)]
       zcur=np.log1p(a);vals=[]
       for k in [i-1,i+1]:
        zn=np.log1p(X[k,j-12:j+13]);gain=np.median(zcur[sh]-zn[sh]);vals.append(np.expm1(zn+gain))
       predtemp=np.mean(vals,axis=0)
       # censored Gaussian peak model: linear baseline + Gaussian peak
       xn=xx/max(np.max(abs(xx)),1e-6);C=(cap-base[12])/max(hgt,1);ys=(obs-base[12])/max(hgt,1);xscale=max(hgt,1);sigN=max(sig/xscale,.01)
       vis=~mask
       def fun(p):
        b0,b1,amp,mu,sd=p;fit=b0+b1*xn+amp*np.exp(-.5*((xn-mu)/sd)**2)
        r=(ys[vis]-fit[vis])/sigN
        return float(.5*np.sum(r*r)-np.sum(log_ndtr((fit[mask]-C)/sigN)))
       best=None
       for initamp in [1,1.5,2.5]:
        fit=minimize(fun,[0,0,initamp,0,.2],method='L-BFGS-B',bounds=[(-2,2),(-2,2),(0,5),(-.5,.5),(.03,1.5)],options={'maxiter':300})
        if best is None or fit.fun<best.fun:best=fit
       b0,b1,amp,mu,sd=best.x;fitnorm=b0+b1*xn+amp*np.exp(-.5*((xn-mu)/sd)**2);predcen=base[12]+xscale*fitnorm
       # use model only at clipped samples; preserve visible samples
       methods={'uncorrected':obs.copy(),'linear':predlin,'neighbor_template':predtemp,'censored_gaussian':np.where(mask,predcen,obs)}
       bl=np.linspace(a[0],a[-1],len(a));true_area=np.trapezoid(np.maximum(a-bl,0),xx);true_height=a.max()-bl[12]
       for name,pred in methods.items():
        e=np.mean(abs(pred[mask]-a[mask]));ph=pred.max()-bl[12];pa=np.trapezoid(np.maximum(pred-bl,0),xx)
        simrows.append({'heat':h,'sequence':i+1,'wavelength_nm':wl[j],'cap_fraction':frac,'method':name,'censored_points':int(mask.sum()),'MAE_censored':float(e),'peak_height_rel_error':float(abs(ph-true_height)/max(abs(true_height),1)),'peak_area_rel_error':float(abs(pa-true_area)/max(abs(true_area),1)),'optimizer_success':bool(best.success) if name=='censored_gaussian' else None})
     print('saturation analyzed',h,'actual_sat',int(sat.sum()),'sim_peaks',len(cand),flush=True)
    pd.DataFrame(summary).to_csv(R/'q1_saturation_structure.csv',index=False);pd.DataFrame(regions).to_csv(R/'q1_saturation_components.csv',index=False);pd.DataFrame(simrows).to_csv(R/'q1_saturation_simulation_allheats.csv',index=False)
    print('\nSTRUCTURE\n',pd.DataFrame(summary).to_string(index=False))
    print('\nSIMULATION\n',pd.DataFrame(simrows).groupby(['cap_fraction','method'])[['MAE_censored','peak_height_rel_error','peak_area_rel_error']].median().round(3).to_string())
    

def main():
    run()

if __name__ == '__main__':
    main()
