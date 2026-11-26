"""Q3 causal cross-heat soft sensing; no gas inputs at publication."""
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from pathlib import Path

GASES = ('CO2', 'CO')
LAGS = np.arange(0, 61, 5)
PENALTIES = (.01, .1, 1.)
PREFIX, WARMUP = 60, 120


def features(wl, x, oxygen, dt):
    """Fixed physical bands, same-frame spectrum, backward differences only."""
    edges = np.linspace(345, 1015, 33)
    bands = np.stack([np.log1p(np.mean(x[:, (wl >= a) & (wl < b)], axis=1))
                      for a, b in zip(edges[:-1], edges[1:])], axis=1)
    if not np.isfinite(bands).all():
        raise ValueError('Invalid spectrum or wavelength support')
    i = np.arange(len(x)); k = np.maximum(0, i-int(10/dt)); j = np.maximum(0, i-int(5/dt))
    q = np.maximum(0, (oxygen-oxygen[k])/np.maximum((i-k)*dt, dt))
    sat = np.mean(x >= 65535, axis=1)
    p = np.column_stack([i*dt/600, (oxygen-oxygen[0])/3000, q/5, (q-q[j])/5])
    b = bands.mean(axis=1)
    s = np.column_stack([bands, bands-bands[j], sat])
    m = i*dt <= PREFIX
    prefix = np.r_[b[m].mean(), b[m].std(), (b-b[j])[m].std(), sat[m].mean(), q[m].mean(), q[m].std(), bands[m].mean(0).reshape(4,8).mean(1)-b[m].mean()]
    contrast = np.stack([np.log1p(np.quantile(x[:, (wl >= a) & (wl < c)], .95, axis=1))-np.log1p(np.median(x[:, (wl >= a) & (wl < c)], axis=1))
                         for a, c in zip(edges[:-1], edges[1:])], axis=1)
    sp = np.column_stack([s, contrast, contrast-contrast[j]])
    return p, s, sp, prefix, b


def prepare(raw, scenario):
    data = {}
    for h, a in raw.items():
        dt = 1. if scenario == 'H1' and h == '6C02057' else .5
        p, s, sp, prefix, brightness = features(a['wl'], a['x'], a['oxygen'], dt)
        gas = a['gas']; idx = np.arange(int(WARMUP/dt), len(gas))
        # Equal detector support for all candidate leads. Backward 5 s differences remove level trend.
        k = int(5/dt); dy = gas[idx]-gas[idx-k]; profile = []
        for d in LAGS:
            j = idx-int(d/dt); dx = brightness[j]-brightness[j-k]
            corr = [np.corrcoef(dx, dy[:,g])[0,1] if dx.std()>1e-12 and dy[:,g].std()>1e-12 else 0. for g in range(2)]
            profile.append(np.abs(corr))
        data[h] = dict(p=p,s=s,sp=sp,prefix=prefix,gas=gas,dt=dt,idx=idx,profile=np.array(profile))
    return data


def lag(data, train, g):
    return int(LAGS[np.argmax(np.mean([data[h]['profile'][:,g] for h in train],axis=0))])


def weighted_stats(xs):
    mu = np.mean([x.mean(0) for x in xs],axis=0)
    sd = np.sqrt(np.mean([((x-mu)**2).mean(0) for x in xs],axis=0))
    return mu, np.where(sd>1e-9,sd,1.)


def normalize(xs, train):
    mu, sd = weighted_stats([xs[h] for h in train])
    return {h:(x-mu)/sd for h,x in xs.items()}


def fold(data, train, test, g, lead):
    hs = list(train)+[test]; pp={}; ss={}; sspeak={}; yy={}; yzero={}; prefix={}
    for h in hs:
        a=data[h]; j=a['idx']-int(lead/a['dt'])
        pp[h]=a['p'][j]; ss[h]=a['s'][j]; sspeak[h]=a['sp'][j]; yy[h]=a['gas'][a['idx'],g]; yzero[h]=a['gas'][j,g]; prefix[h]=a['prefix']
    p=normalize(pp,train); s=normalize(ss,train)
    cov=np.mean([s[h].T@s[h]/len(s[h]) for h in train],axis=0)
    _,v=np.linalg.eigh(cov); z={h:s[h]@v[:,-8:] for h in hs}
    sp=normalize(sspeak,train); cp=np.mean([sp[h].T@sp[h]/len(sp[h]) for h in train],axis=0)
    _,vp=np.linalg.eigh(cp); zp={h:sp[h]@vp[:,-8:] for h in hs}
    # Identical quadratic process basis in nonlinear process-only and state interaction models.
    p2={h:np.column_stack([p[h],*[p[h][:,i]*p[h][:,j] for i in range(4) for j in range(i,4)]]) for h in hs}
    designs={'process':p,'process_nonlinear':p2,'spectral':z,
             'fusion':{h:np.c_[p2[h],z[h]] for h in hs},
             'spectral_peak':zp,'fusion_peak':{h:np.c_[p2[h],zp[h]] for h in hs},
             'state':{h:np.c_[p2[h],z[h],z[h]*p[h][:,1:2],z[h]*p[h][:,2:3]] for h in hs}}
    designs={name:normalize(x,train) for name,x in designs.items()}
    designs={name:{h:np.c_[np.ones(len(x[h])),x[h]] for h in hs} for name,x in designs.items()}
    u=np.stack([prefix[h] for h in train]); mu=u.mean(0); sd=u.std(0); sd=np.where(sd>1e-9,sd,1.)
    u=(u-mu)/sd; ut=(prefix[test]-mu)/sd
    return designs, s, yy, yzero, u, ut


def ridge_parts(x,y,train):
    return np.stack([x[h].T@x[h]/len(y[h]) for h in train]), np.stack([x[h].T@y[h]/len(y[h]) for h in train])


def solve(a,b,pen):
    reg=np.eye(len(b))*pen; reg[0,0]=0
    return np.linalg.solve(a+reg,b)


def pls(s,y,train,test,r):
    # PLS1 deflation in a heat-equal weighted inner product.
    mu=np.mean([y[h].mean() for h in train]); X=np.concatenate([s[h]/np.sqrt(len(y[h])*len(train)) for h in train]); Y=np.concatenate([(y[h]-mu)/np.sqrt(len(y[h])*len(train)) for h in train]); W=[]; P=[]; Q=[]
    for _ in range(r):
        w=X.T@Y; norm=np.linalg.norm(w)
        if norm<1e-12:break
        w/=norm; t=X@w; den=t@t
        if den<1e-12:break
        p=X.T@t/den; q=Y@t/den; X-=np.outer(t,p); Y-=t*q; W.append(w);P.append(p);Q.append(q)
    if not W:return np.full(len(s[test]),mu)
    w=np.array(W).T; p=np.array(P).T
    return mu+s[test]@w@np.linalg.pinv(p.T@w)@np.array(Q)


def predict_candidates(data,train,test,g,lead):
    d,s,y,y0,u,ut=fold(data,train,test,g,lead); n=len(train); pred={}
    pred[('mean',0)]=np.full(len(y[test]),np.mean([y[h].mean() for h in train]))
    for r in (4,8):pred[('spectral_pls',r)]=pls(s,y,train,test,r)
    dist=np.mean((u-ut)**2,axis=1); sims=np.exp(-(dist-dist.min())/2); sims/=sims.sum()
    for name,x in d.items():
        A,B=ridge_parts(x,y,train); a=A.mean(0); b=B.mean(0)
        for pen in PENALTIES:
            coef=solve(a,b,pen); base=x[test]@coef; pred[(name,pen)]=base
            if name in ('process','process_nonlinear','fusion_peak'):
                own=np.mean([x[test]@solve(A[i],B[i],pen) for i in range(n)],axis=0)
                pred[(name+'_ensemble',pen)]=.5*base+.5*own
            if name!='fusion':continue
            _,B0=ridge_parts(x,y0,train); pred[('zero_alignment',pen)]=x[test]@solve(a,B0.mean(0),pen)
            # Historical furnace models and affine calibration from cross-fitted global predictions.
            individual=[]; ab=[]; scale=max(1.,np.std(np.concatenate([y[h] for h in train])))
            for i,h in enumerate(train):
                own=solve(A[i],B[i],pen); individual.append(x[test]@own)
                cf=x[h]@solve((A.sum(0)-A[i])/(n-1),(B.sum(0)-B[i])/(n-1),pen)
                C=np.c_[np.ones(len(cf)),cf/scale]; target=(y[h]-cf)/scale
                delta=np.linalg.solve(C.T@C/len(cf)+np.diag([.1,.1]),C.T@target/len(cf)); ab.append([delta[0]*scale,np.clip(1+delta[1],.2,3.)])
            pred[('individual_ensemble',pen)]=.5*base+.5*np.mean(individual,axis=0)
            alpha,beta=sims@np.array(ab)
            for shrink in (.5,1.):pred[('alpha_beta',(pen,shrink))]=(1-shrink)*base+shrink*(alpha+beta*base)
            for K in (2,3):
                km=KMeans(n_clusters=K,n_init=10,random_state=42).fit(u); labels=km.labels_; local=[]; support=[]; distance=[]
                for k in range(K):
                    ids=np.flatnonzero(labels==k); m=len(ids)
                    if not m:continue
                    local.append(x[test]@solve(A[ids].mean(0),B[ids].mean(0),pen)); support.append(m/(m+2)); distance.append(np.mean((ut-km.cluster_centers_[k])**2))
                local=np.array(local); support=np.array(support); distance=np.array(distance)
                k=np.argmin(distance); pred[(f'hard_moe{K}',pen)]=(1-support[k])*base+support[k]*local[k]
                for temp in (.5,1.,2.):
                    weight=support*np.exp(-(distance-distance.min())/(2*temp**2)); weight/=weight.sum()
                    strength=.75*np.exp(-distance.min()/8)*float(weight@support)
                    pred[(f'soft_moe{K}',(pen,temp))]=(1-strength)*base+strength*(weight@local)
    return pred,y[test]


def run(raw, scenarios=('H0','H1')):
    scores=[]; predictions=[]; lagrows=[]; choices=[]
    for scenario in scenarios:
        data=prepare(raw,scenario); hs=sorted(data)
        for g,gas in enumerate(GASES):
            for test in hs:
                train=[h for h in hs if h!=test]; inner={}
                for val in train:
                    tr=[h for h in train if h!=val]; lead=lag(data,tr,g)
                    ps,y=predict_candidates(data,tr,val,g,lead)
                    for key,v in ps.items():inner.setdefault(key,[]).append(np.mean(np.abs(np.maximum(0,v)-y)))
                best={}
                for key,err in inner.items():
                    name=key[0]; e=np.mean(err)
                    if name not in best or e<best[name][0]:best[name]=(e,key)
                lead=lag(data,train,g); ps,y=predict_candidates(data,train,test,g,lead); a=data[test]; det=a['idx']*a['dt']
                lagrows.append([scenario,gas,test,lead,lead in (0,60)])
                for name,(inner_mae,key) in best.items():
                    v=ps[key]; vh=np.maximum(0,v); err=vh-y; den=np.sum((y-y.mean())**2)
                    scores.append([scenario,gas,test,name,lead,len(y),np.abs(err).mean(),np.sqrt(np.mean(err**2)),1-np.sum(err**2)/den if den>0 else np.nan,np.mean(v<0),inner_mae,str(key[1])])
                    choices.append([scenario,gas,test,name,str(key[1]),inner_mae])
                    predictions.extend(zip([scenario]*len(y),[gas]*len(y),[test]*len(y),[name]*len(y),det-lead,det,y,v,vh))
                print(scenario,gas,test,lead,flush=True)
    scores=pd.DataFrame(scores,columns=['scenario','gas','heat','model','lead_s','n','mae','rmse','r2','negative_share','inner_mae','parameter'])
    pred=pd.DataFrame(predictions,columns=['scenario','gas','heat','model','publish_s','detector_s','observed','raw_prediction','prediction'])
    macro=scores.groupby(['scenario','gas','model']).agg(heats=('heat','nunique'),mae=('mae','mean'),rmse=('rmse','mean'),r2=('r2','mean'),worst_mae=('mae','max')).reset_index()
    return scores,pred,macro,pd.DataFrame(lagrows,columns=['scenario','gas','heat','lead_s','boundary']),pd.DataFrame(choices,columns=['scenario','gas','heat','model','parameter','inner_mae'])


def save(raw,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    scores,pred,macro,lags,choices=run(raw)
    for name,df in [('fold_metrics',scores),('macro_metrics',macro),('lags',lags),('inner_choices',choices)]:df.to_csv(out/(name+'.csv'),index=False,float_format='%.8g')
    pred.to_csv(out/'predictions.csv.gz',index=False,compression={'method':'gzip','mtime':0},float_format='%.8g')
    print(macro.to_string(index=False),flush=True)
