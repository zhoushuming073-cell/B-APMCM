"""Equal-support correlation and exhaustive remote circular-shift diagnostic."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from .delay import series
from .reconstruct import TAB,GASES

def full_profiles(curves):
    stats=[];profiles=[]
    for (scen,heat),frame in curves.groupby(['scenario','heat']):
        t,q,_,yy=series(frame,8)
        x=np.gradient(gaussian_filter1d(q,1.5))[:-30]
        x=x-x.mean()
        m=len(x)
        for g,gas in enumerate(GASES):
            y=np.gradient(gaussian_filter1d(yy[:,g],1.5))
            n=len(y)
            windows=np.lib.stride_tricks.sliding_window_view(np.r_[y,y],m)[:n]
            centered=windows-windows.mean(axis=1)[:,None]
            denom=np.linalg.norm(centered,axis=1)*np.linalg.norm(x)
            r=np.divide(centered@x,denom,out=np.zeros(n),where=denom>1.e-10)
            lags=np.arange(31)
            k=int(np.argmax(abs(r[lags])))
            shifts=np.arange(45,n-44)
            null=abs(r[(shifts[:,None]+lags)%n]).max(axis=1)
            stats.append(dict(scenario=scen,heat=heat,gas=gas,lag_s=2*k,r=r[k],boundary=k in (0,30),
                              p=(1+(null>=abs(r[k])).sum())/(len(null)+1),n_shifts=len(shifts)))
            profiles.extend((scen,heat,gas,2*l,r[l]) for l in lags)
    df=pd.DataFrame(stats)
    for scen,g in df.groupby('scenario'):
        p=g.p.to_numpy()
        order=p.argsort()
        adj=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
        df.loc[g.index[order],'bh_q']=np.minimum(adj,1)
    df.to_csv(TAB/'global_lags.csv',index=False)
    pd.DataFrame(profiles,columns=['scenario','heat','gas','lag_s','r']).to_csv(TAB/'global_profiles.csv',index=False)
