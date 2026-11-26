"""Offline continuous reconstruction, internal-gap validation, and dose associations."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator, PPoly
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT/'tables'
GASES = ['CO2','CO']
CHANNELS = ['O',*GASES]


def flags(t,v,var,dt):
    # A brief reversal must disagree with both flanks, which themselves must agree.
    flag = np.zeros(len(v),bool)
    if var=='O':
        # A counter cannot decrease in these no-reset records. Of the two rows
        # bracketing a decrease, flag the larger two-sided interpolation residual.
        for j in np.flatnonzero(np.diff(v)<0):
            cand=np.array([i for i in (j,j+1) if 0<i<len(v)-1 and t[i+1]-t[i-1]<=2*dt+1.e-8],dtype=int)
            if len(cand):
                residual=abs(v[cand]-(v[cand-1]+v[cand+1])/2)
                flag[cand[np.argmax(residual)]]=True
        return flag
    noise_floor = .2
    from numpy.lib.stride_tricks import sliding_window_view
    w=sliding_window_view(v,9)
    flanks=w[:,[0,1,2,3,5,6,7,8]]
    left,right=np.median(flanks[:,:4],axis=1),np.median(flanks[:,4:],axis=1)
    med=np.median(flanks,axis=1)
    limit=np.maximum(noise_floor,6*1.4826*np.median(abs(flanks-med[:,None]),axis=1))
    center=w[:,4]
    distance=np.minimum.reduce([abs(center-w[:,3]),abs(center-w[:,5]),abs(center-left),abs(center-right)])
    flag[4:-4]=((t[8:]-t[:-8]<=8*dt+1.e-8)&(abs(left-right)<limit)&(abs(w[:,3]-w[:,5])<limit)
                &((center-w[:,3])*(center-w[:,5])>0)&(distance>limit))
    return flag


class AveragedCounter:
    """Box-average a monotone polygon, reflect it about its endpoints to preserve both counter readings."""
    def __init__(self,t,v,width):
        self.h=width/2
        self.t=t
        self.v=np.clip(np.maximum.accumulate(v),v[0],v[-1])
        self.linear=PPoly(np.array([np.diff(self.v)/np.diff(t),self.v[:-1]]),t)
        self.integral=self.linear.antiderivative()
        start,end=self.average(np.array([t[0],t[-1]]))
        self.start=start
        self.scale=(v[-1]-v[0])/(end-start) if end>start else 0.

    def primitive(self,x):
        lo,hi=self.t[0],self.t[-1]
        reflected=np.where(x<lo,2*lo-x,np.where(x>hi,2*hi-x,x))
        return self.integral(reflected)+np.minimum(x-lo,0)*2*self.v[0]+np.maximum(x-hi,0)*2*self.v[-1]

    def average(self,x):
        return (self.primitive(x+self.h)-self.primitive(x-self.h))/(2*self.h)

    def __call__(self,x):
        return self.v[0]+self.scale*(self.average(np.asarray(x))-self.start)

    def derivative(self):
        def flow(x):
            x=np.asarray(x)
            def extended(z):
                lo,hi=self.t[0],self.t[-1]
                mirrored=np.where(z<lo,2*lo-z,np.where(z>hi,2*hi-z,z))
                v=np.interp(mirrored,self.t,self.v)
                return np.where(z<lo,2*self.v[0]-v,np.where(z>hi,2*self.v[-1]-v,v))
            return self.scale*(extended(x+self.h)-extended(x-self.h))/(2*self.h)
        return flow


def build(t,v,var,method,dt):
    if method != 'linear':
        keep = ~flags(t,v,var,dt)
        t,v=t[keep],v[keep]
    if method.startswith('a'):
        return AveragedCounter(t,v,float(method[1:]))
    if method in ('linear','clean_linear'):
        return PPoly(np.array([np.diff(v)/np.diff(t),v[:-1]]),t,extrapolate=False)
    if method != 'pchip':
        width = float(method[1:])
        knots=np.unique(np.r_[np.arange(t[0],t[-1],width),t[-1]])
        # Only observed samples contribute. Knots falling inside gaps are skipped.
        centers,values=[t[0]],[v[0]]
        for k in knots[1:-1]:
            sel=abs(t-k)<=width/2
            if sel.sum()>=2:
                x,y=t[sel],v[sel]
                centers.append(x.mean())
                values.append(y.mean())
        t,v=np.r_[centers,t[-1]],np.r_[values,v[-1]]
    if var=='O':
        v=np.minimum(v[-1],np.maximum.accumulate(v))
    return PchipInterpolator(t,np.maximum(v,0),extrapolate=False)


def metrics(y,z):
    e=z-y
    return dict(mae=np.mean(abs(e)),rmse=np.sqrt(np.mean(e**2)),
                r2=1-np.sum(e**2)/np.sum((y-y.mean())**2) if np.ptp(y)>0 else np.nan)


def candidates(var):
    return ['linear','pchip','k4','k8','k16','a4','a8','a16'] if var=='O' else ['linear','clean_linear','pchip','k2','k4','k8']


def reconstruct(data):
    scores=[]; curves=[]; audit=[]; associations=[]; replacements=[]; params=[]
    for scen in ['H0','H1']:
        hdata={}
        for heat,g in data.groupby('heat'):
            a=int(g.attachment.iloc[0])
            dt=1. if scen=='H1' and a==5 else .5
            t=np.arange(len(g))*dt
            hdata[heat]=(a,dt,t,g)
            # Same row blocks in H0/H1. 3 interleaved sets, 12 internal blocks total.
            for var in CHANNELS:
                v=g[var].to_numpy()
                for fold in range(3):
                    blocks=np.rint(np.linspace(.06,.92,12)[fold::3]*(len(t)-1)).astype(int)
                    hide=np.zeros(len(t),bool)
                    for start,n in zip(blocks,[8,16,8,16]):
                        hide[start:start+n]=True
                    for method in candidates(var):
                        f=build(t[~hide],v[~hide],var,method,dt)
                        scores.append(dict(scenario=scen,heat=heat,channel=var,fold=fold+1,method=method,
                                           n=hide.sum(),sum_y=v[hide].sum(),sum_y2=(v[hide]**2).sum(),**metrics(v[hide],f(t[hide]))))
        cv=pd.DataFrame(scores).query('scenario == @scen')
        for heat,(a,dt,t,g) in hdata.items():
            out=dict(scenario=scen,heat=heat,t=t,row=g.row.to_numpy(),attachment=a)
            for var in CHANNELS:
                v=g[var].to_numpy()
                # Continuous flow requires the derivative to be continuous, so exclude
                # the polygonal O and the unregularized staircase interpolant from deployment.
                eligible=['k4','k8','k16','a4','a8','a16'] if var=='O' else [m for m in candidates(var) if m!='linear']
                pool=cv[(cv.heat!=heat)&(cv.channel==var)&cv.method.isin(eligible)]
                mean=pool.groupby('method').mae.mean()
                method=mean.idxmin()
                params.append(dict(scenario=scen,heat=heat,channel=var,method=method))
                f=build(t,v,var,method,dt)
                out[var]=f(t)
                out[var+'_raw']=v
                marked=np.where(flags(t,v,var,dt))[0]
                for i in marked:
                    replacements.append(dict(scenario=scen,heat=heat,channel=var,row=i+1,t=t[i],raw=v[i],
                                             reconstructed=f(t[i]),used=method!='linear'))
                if var=='O':
                    out['q']=f.derivative()(t)
                else:
                    F=f.antiderivative()
                    out['E_'+var]=F(t)-F(0.)
                audit.append(dict(scenario=scen,heat=heat,channel=var,method=method,n=len(v),
                                  n_flag=len(marked),applied_flags=len(marked) if method!='linear' else 0,
                                  start=f(t[0]),end=f(t[-1]),duration_s=t[-1],**metrics(v,f(t))))
            frame=pd.DataFrame(out)
            curves.append(frame)
            o=frame.O.to_numpy()
            for var in ['CO2','CO','E_CO2','E_CO']:
                y=frame[var].to_numpy()
                ix=np.arange(0,len(t),round(30/dt))
                # Disjoint blocks: oxygen increments vs gas averages or emission increments.
                block=np.diff(y[ix]) if var.startswith('E_') else np.array([np.mean(y[l:r]) for l,r in zip(ix,ix[1:])])
                dx=np.diff(o[ix])
                associations.append(dict(scenario=scen,heat=heat,outcome=var,n_blocks=len(dx),
                    rho_level=spearmanr(o,y).statistic,rho_blocks=spearmanr(dx,block).statistic,
                    end=y[-1],oxygen_net=o[-1]-o[0]))
    tables={'reconstruction_cv':scores,'reconstruction_choice':params,'curve_audit':audit,'associations':associations}
    for name,rows in tables.items():
        pd.DataFrame(rows).to_csv(TAB/(name+'.csv'),index=False)
    pd.DataFrame(replacements,columns=['scenario','heat','channel','row','t','raw','reconstructed','used']).to_csv(TAB/'spike_flags.csv',index=False)
    all_curves=pd.concat(curves,ignore_index=True)
    all_curves.to_csv(TAB/'curves.csv.gz',index=False,float_format='%.10g',compression={'method':'gzip','mtime':0})
    return all_curves
