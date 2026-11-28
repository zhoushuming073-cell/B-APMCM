import logging
from scipy.ndimage import median_filter, gaussian_filter, maximum_filter1d
from scipy.interpolate import RectBivariateSpline
from apmcm.common import *

# The two residual tests reject an isolated point only when both axes disagree.
def repair(x):
    x=x.copy()
    missing=~np.isfinite(x)
    if missing.any():
        for j in np.flatnonzero(missing.any(axis=0)):
            ok=~missing[:,j]
            if not ok.any(): raise ValueError("Entire wavelength missing")
            x[:,j]=np.interp(np.arange(len(x)),np.flatnonzero(ok),x[ok,j])
    time=median_filter(x,size=(5,1),mode="nearest")
    wave=median_filter(x,size=(1,5),mode="nearest")
    rt=x-time; rw=x-wave
    scales=[]
    for resid in [rt,rw]:
        local=median_filter(abs(resid),size=(3,5),mode="nearest")
        floor=np.maximum(1,np.median(abs(resid),axis=0))
        scales.append(1.4826*np.maximum(local,floor))
    bad=(np.abs(rt)>8*scales[0]) & (np.abs(rw)>8*scales[1]) & (rt*rw>0)
    x[bad]=((time+wave)/2)[bad]
    return x,bad,missing

OPTIONS=[(1,0.),(3,0.),(1,.55),(3,.55),(3,1.)]
def surface(t,w,x,k,sigma):
    z=gaussian_filter(x,sigma,mode="nearest") if sigma else x
    return RectBivariateSpline(t,w,z,kx=k,ky=k,s=0)

def run(heats):
    cv=[]; metrics=[]
    upper=max(float(np.nanpercentile(h["x"],99.85)) for h in heats.values())
    for ix,(heat,h) in enumerate(heats.items()):
        x,t,w=h["x"],h["t"],h["wl"]
        it=np.arange(len(t));iw=np.arange(len(w))
        kt=(it%8!=3);kw=(iw%8!=3)
        train=repair(x[np.ix_(kt,kw)])[0]
        # Split by disjoint time blocks; the withheld grid lines are absent during every fit.
        rng=np.random.default_rng(SEED+ix)
        held_t=it[(~kt)&(it>4)&(it<len(t)-5)]
        valid_t,test_t=np.array_split(rng.permutation(held_t),2)
        sets=[]
        for rows in (valid_t,test_t):
            a=rng.choice(rows,14000); b=rng.integers(4,len(w)-4,14000)
            # Also score withheld wavelength lines on retained time slices (disjoint by half).
            rt=rng.choice(it[kt & ((it%16<8) if rows is valid_t else (it%16>=8))],14000)
            rb=rng.choice(iw[~kw],14000)
            a=np.r_[a,rt];b=np.r_[b,rb]
            ok=np.isfinite(x[a,b]);sets.append((a[ok],b[ok]))
        local=[]
        for k,sigma in OPTIONS:
            f=surface(t[kt],w[kw],train,k,sigma)
            for split,(a,b) in zip(["validation","test"],sets):
                err=f.ev(t[a],w[b])-x[a,b]
                row=dict(heat=heat,degree=k,sigma=sigma,split=split,n=len(err),
                         mae=np.mean(abs(err)),rmse=np.sqrt(np.mean(err**2)))
                cv.append(row);local.append(row)
        choice=min((r for r in local if r["split"]=="validation"),key=lambda r:r["mae"])
        k,sigma=choice["degree"],choice["sigma"]
        clean,bad,missing=repair(x)
        fitted=gaussian_filter(clean,sigma,mode="nearest") if sigma else clean.copy()
        np.savez(MODEL/f"spectrum_{heat}.npz",clean=clean,nodes=fitted,flag=bad,missing=missing,degree=k)
        f=surface(t,w,fitted,k,0)
        peaks=(clean==maximum_filter1d(clean,9,axis=1)) & (clean>np.quantile(clean,.95))
        ratio=np.median(fitted[peaks]/np.maximum(clean[peaks],1)) if peaks.any() else np.nan
        rough=lambda z:float(np.mean(abs(np.diff(z,n=2,axis=1))))
        result=next(r for r in local if r["split"]=="test" and (r["degree"],r["sigma"])==(k,sigma))
        metrics.append(dict(heat=heat,cells=x.size,outliers=bad.sum(),fraction=bad.mean(),missing=missing.sum(),
                            degree=k,sigma=sigma,test_mae=result["mae"],test_rmse=result["rmse"],
                            peak_ratio=ratio,rough_raw=rough(x),rough_fit=rough(fitted)))
        fig=plt.figure(figsize=(8,5));ax=fig.add_subplot(projection="3d")
        ts=np.linspace(t[0],t[-1],100);ws=np.linspace(w[0],w[-1],145)
        tt,ww=np.meshgrid(ts,ws,indexing="ij");zz=f(ts,ws)
        surf=ax.plot_surface(tt,ww,zz,cmap="viridis",vmin=0,vmax=upper,linewidth=0,rcount=100,ccount=145,rasterized=True)
        ax.set(xlabel="Time (s)",ylabel="Wavelength (nm)",zlabel="Intensity (W/cm²)",title=f"Reconstructed flame spectrum | {heat}")
        fig.colorbar(surf,ax=ax,shrink=.58,label="Intensity (W/cm²)");figure(fig,"Q1_surface_"+heat)
        delta=np.abs(x-clean)
        a,b=np.unravel_index(np.nanargmax(delta),delta.shape)
        sl=slice(max(0,b-25),min(len(w),b+26))
        fig,axes=plt.subplots(1,2,figsize=(9,3.1))
        axes[0].plot(w[sl],x[a,sl],label="Raw");axes[0].plot(w[sl],fitted[a,sl],label="Reconstruction")
        axes[0].set(xlabel="Wavelength (nm)",ylabel="Intensity (W/cm²)",title=f"At {t[a]:.1f} s")
        sl=slice(max(0,a-25),min(len(t),a+26))
        axes[1].plot(t[sl],x[sl,b],label="Raw");axes[1].plot(t[sl],fitted[sl,b],label="Reconstruction")
        axes[1].set(xlabel="Time (s)",title=f"At {w[b]:.3f} nm")
        for ax in axes: ax.legend()
        fig.suptitle(f"Largest repair | {heat}");fig.tight_layout();figure(fig,"Q1_repair_"+heat)
        logging.info("Spectrum %s: %d flags, degree %s sigma %s",heat,bad.sum(),k,sigma)
    table("q1_validation",cv);table("q1_metrics",metrics)
