"""Compare delay laws on identical held-out gas readings, with nested heat splits."""
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from reconstruct import TAB, GASES

QREF=5.
LAGS=[0,6,12,18,24,30,36,42,48,54,60]
BETAS=[0.,np.log(1.5),np.log(2.),-np.log(1.5),-np.log(2.)]
CONFIG=[('stage',0.,0.)]+[('shift',float(d),b) for d in LAGS for b in (BETAS if d else [0.])]


def series(frame,window):
    t0=frame.t.to_numpy()
    t=np.arange(0,t0[-1]+1.e-8,2.)
    oxygen=frame.O_raw.to_numpy()
    at=lambda x: np.interp(x,t0,oxygen)
    lo=np.maximum(0,t-window/2)
    hi=np.minimum(t0[-1],t+window/2)
    q=(at(hi)-at(lo))/(hi-lo)
    u=(at(t)-oxygen[0])/3000
    stage=np.column_stack([np.ones(len(t)),t/600,u,*[np.maximum(u-k,0) for k in (.15,.4,.7,.9)]])
    y=np.column_stack([np.interp(t,t0,frame[g+'_raw']) for g in GASES])
    return t,q,stage,y


def design(series,config):
    t,q,stage,_=series
    kind,d,b=config
    valid=t>=120
    if kind=='stage':
        return stage[valid],None
    speed=np.exp(b*np.clip((q-QREF)/QREF,-1,1))
    clock=cumulative_trapezoid(speed,t,initial=0)
    need=clock[valid]-d
    assert need.min()>=-1.e-7
    source=np.interp(need,clock,t)
    aligned=np.interp(source,t,q)/QREF
    return np.column_stack([stage[valid],aligned,aligned**2]),t[valid]-source


def eligible(name):
    if name=='stage':
        return [0]
    if name=='zero':
        return [1]
    if name=='fixed':
        return [i for i,(k,d,b) in enumerate(CONFIG) if k=='shift' and b==0]
    sign=1 if name=='decreasing' else -1
    return [i for i,(k,d,b) in enumerate(CONFIG) if k=='shift' and b*sign>=0]


def delay_study(curves):
    scores=[];profiles=[];bands=[];chosen=[];preds=[]
    for scen in ['H0','H1','H0_without_5']:
        subset=curves[curves.scenario==('H0' if scen=='H0_without_5' else scen)]
        if scen=='H0_without_5':
            subset=subset[subset.attachment!=5]
        group=list(subset.groupby('heat'))
        heats=[h for h,_ in group]
        for window in [4,8,16] if scen=='H0' else [8]:
            series_set=[series(f,window) for _,f in group]
            n=len(group)
            # for a specified configuration: loss of validation heat v with outer heat h excluded
            inner=np.full((len(CONFIG),n,n,2),np.nan)
            outer=np.zeros((len(CONFIG),n,2))
            saved={}
            for c,config in enumerate(CONFIG):
                xys=[]
                for s in series_set:
                    x,_=design(s,config)
                    y=s[3][s[0]>=120]
                    xys.append((x,y))
                gram=np.array([x.T@x/len(x) for x,y in xys])
                rhs=np.array([x.T@y/len(x) for x,y in xys])
                penalty=np.eye(gram.shape[-1])*.001
                penalty[0,0]=0
                A,B=gram.sum(axis=0),rhs.sum(axis=0)
                for h in range(n):
                    x,y=xys[h]
                    coef=np.linalg.solve(A-gram[h]+(n-1)*penalty,B-rhs[h])
                    yhat=x@coef
                    saved[c,h]=yhat
                    outer[c,h]=np.mean(abs(yhat-y),axis=0)
                    for v in range(h):
                        coef=np.linalg.solve(A-gram[h]-gram[v]+(n-2)*penalty,B-rhs[h]-rhs[v])
                        for a,j in [(h,v),(v,h)]:
                            xx,yy=xys[j]
                            inner[c,a,j]=np.mean(abs(xx@coef-yy),axis=0)
                for g,gas in enumerate(GASES):
                    profiles.append(dict(scenario=scen,window=window,gas=gas,kind=config[0],d=config[1],
                                         beta=config[2],macro_mae=outer[c,:,g].mean()))
            avg=np.nanmean(inner,axis=2)
            for model in ['stage','zero','fixed','decreasing','reverse']:
                ids=eligible(model)
                for g,gas in enumerate(GASES):
                    # Full-data calibration using leave-one-heat-out, never used for reporting outer accuracy.
                    best=ids[np.argmin(outer[ids,:,g].mean(axis=1))]
                    _,d,b=CONFIG[best]
                    chosen.append(dict(scenario=scen,window=window,gas=gas,model=model,d=d,beta=b))
                    if window==8 and model in ['fixed','decreasing']:
                        qq=np.concatenate([s[1][s[0]>=120] for s in series_set])
                        for ref in [3.,4.,5.,6.,8.,10.]:
                            count=sum(np.any(abs(s[1][s[0]>=120]-ref)<.25) for s in series_set)
                            bands.append(dict(scenario=scen,model=model,gas=gas,q=ref,d=d,beta=b,
                                              local_constant_tau=d/np.exp(b*np.clip((ref-QREF)/QREF,-1,1)),
                                              observed_points=np.sum(abs(qq-ref)<.25),heats=count))
                    for h,heat in enumerate(heats):
                        c=ids[np.argmin(avg[ids,h,g])]
                        _,d,b=CONFIG[c]
                        s=series_set[h]
                        y=s[3][s[0]>=120,g]
                        yh=saved[c,h][:,g]
                        scores.append(dict(scenario=scen,window=window,heat=heat,gas=gas,model=model,
                                           d=d,beta=b,n=len(y),mae=outer[c,h,g],
                                           rmse=np.sqrt(np.mean((y-yh)**2))))
                        if window==8:
                            valid=s[0]>=120
                            for ti,obs,est in zip(s[0][valid],y,yh):
                                preds.append((scen,heat,gas,model,ti,obs,est))
    for name,rows in [('delay_folds',scores),('delay_profiles',profiles),('delay_calibration',chosen),('delay_by_flow',bands)]:
        pd.DataFrame(rows).to_csv(TAB/(name+'.csv'),index=False)
    pd.DataFrame(preds,columns=['scenario','heat','gas','model','t','raw','estimated']).to_csv(
        TAB/'delay_predictions.csv.gz',index=False,float_format='%.8g',compression={'method':'gzip','mtime':0})


def local_profiles(curves):
    rows=[]
    for (scen,heat),frame in curves.groupby(['scenario','heat']):
        t,q,_,yy=series(frame,8)
        x=np.gradient(gaussian_filter1d(q,1.5))
        for g,gas in enumerate(GASES):
            y=np.gradient(gaussian_filter1d(yy[:,g],1.5))
            for j in range(0,len(t)-44,45):
                a=x[j:j+45][:30]
                b=y[j:j+45]
                if a.std()<1.e-7 or b.std()<1.e-7: continue
                cc=[]
                for lag in range(16):
                    v=b[lag:lag+30]
                    cc.append(np.corrcoef(a,v)[0,1] if v.std()>1.e-7 else 0.)
                k=np.argmax(abs(np.array(cc)))
                rows.append(dict(scenario=scen,heat=heat,gas=gas,t=t[j],q=q[j:j+30].mean(),lag=2*k,
                                 r=cc[k],boundary=k in (0,15)))
    df=pd.DataFrame(rows)
    df.to_csv(TAB/'local_windows.csv',index=False)
    outs=[]
    for (scen,gas),h in df.groupby(['scenario','gas']):
        bounds=h.q.quantile([1/3,2/3]).to_numpy()
        for label,grp in h.groupby(pd.cut(h.q,[-np.inf,*bounds,np.inf],labels=['low','middle','high']),observed=True):
            outs.append(dict(scenario=scen,gas=gas,band=label,n_windows=len(grp),n_heats=grp.heat.nunique(),
                             q_min=grp.q.min(),q_max=grp.q.max(),q_median=grp.q.median(),lag_median=grp.lag.median(),
                             lag_q25=grp.lag.quantile(.25),lag_q75=grp.lag.quantile(.75),boundary_rate=grp.boundary.mean()))
    pd.DataFrame(outs).to_csv(TAB/'empirical_flow_bins.csv',index=False)
