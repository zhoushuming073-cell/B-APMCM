"""Reproduce native-grid screening and conservative Q1 template scenarios."""
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d as uf, maximum_filter, median_filter

from apmcm.common.paths import DT, PROC

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


def reconstruct(out, heats=()):
    scenarios = out / "scenarios"
    tables = out / "tables"
    scenarios.mkdir(parents=True, exist_ok=True)
    tables.mkdir(exist_ok=True)
    files = {p.stem[4:]: p for p in sorted(PROC.glob("raw_*.npz"))}
    selected = sorted(set(heats) or files)
    if not selected or not set(selected) <= files.keys():
        raise ValueError("Requested spectra are missing from the official extraction cache")
    rows = []
    for heat in selected:
        with np.load(files[heat]) as source:
            x = source["x"].astype(np.float32)
            wl, seq = source["wl"], source["sequence"]
        saturated = x >= 65535
        train = x[:int(.55 * len(x))]
        scales = fit_scales(train)
        calibration = score(train, scales)
        valid = np.ones(train.shape, bool)
        valid[:15] = valid[-15:] = False
        valid[:, :40] = valid[:, -40:] = False
        threshold = float(np.quantile(calibration[valid], .999))
        full = score(x, scales)
        valid = np.ones(x.shape, bool)
        valid[:15] = valid[-15:] = False
        valid[:, :40] = valid[:, -40:] = False
        screen = (full > threshold) & valid & ~saturated
        value, difference, agree, gap = options(x, screen)
        approved = screen & agree & (difference > .12) & (difference > 2 * gap) & np.isfinite(value)
        scenario = x.copy()
        scenario[approved] = np.expm1(value[approved])
        scenario, upper_estimates = repair_sat(scenario, saturated)
        np.savez_compressed(scenarios / f"{heat}_scenario.npz", raw=x, scenario=scenario,
                            wl=wl, seq=seq, score=full, screen=screen, approved=approved,
                            saturated=saturated, saturation_repaired=upper_estimates)
        n_screen, n_approved = int(screen.sum()), int(approved.sum())
        n_upper, n_est = int(saturated.sum()), int(upper_estimates.sum())
        rows.append(dict(heat=heat, n_time=len(x), screen_candidates=n_screen,
                         screen_rate=n_screen / x.size, approved_anomaly_repairs=n_approved,
                         approved_rate=n_approved / x.size, saturation_points=n_upper,
                         saturation_repaired=n_est, saturation_unresolved=n_upper - n_est,
                         total_changed=n_approved + n_est, threshold=threshold))
        i, j = np.where(screen)
        ix = np.argsort(full[i, j])[::-1][:500]
        i, j = i[ix], j[ix]
        pd.DataFrame(dict(sequence=seq[i], time_s=(seq[i] - seq[0]) * DT, wavelength_nm=wl[j],
                          raw=x[i, j], score=full[i, j], approved_repair=approved[i, j],
                          scenario=scenario[i, j], saturated=saturated[i, j])).to_csv(
            tables / f"{heat}_candidate_top500.csv", index=False)
        print(f"{heat}: {n_screen} candidates, {n_approved} approved", flush=True)
    pd.DataFrame(rows).to_csv(tables / "q1_full_reconstruction_summary.csv", index=False)
