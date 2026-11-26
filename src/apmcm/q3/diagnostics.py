"""Audit equal targets; compute paired heat-level effects and diagnostic charts."""
from pathlib import Path
import numpy as np
import pandas as pd
from .model import features, LAGS, WARMUP


def lag_control(raw, out):
    rows=[]; rng=np.random.default_rng(20261127)
    for sc in ('H0','H1'):
        for gi, gas in enumerate(('CO2','CO')):
            xs=[]; ys=[]; dts=[]
            for h,a in sorted(raw.items()):
                dt=1. if sc=='H1' and h=='6C02057' else .5
                _,_,_,_,b=features(a['wl'],a['x'],a['oxygen'],dt)
                i=np.arange(int(WARMUP/dt),len(b)); k=int(5/dt)
                v=np.stack([b[i-int(d/dt)]-b[i-int(d/dt)-k] for d in LAGS]); y=a['gas'][i,gi]-a['gas'][i-k,gi]
                v=v-v.mean(1,keepdims=True); v=v/np.maximum(np.linalg.norm(v,axis=1,keepdims=True),1e-12)
                y=y-y.mean(); y=y/max(np.linalg.norm(y),1e-12)
                xs.append(v);ys.append(y);dts.append(dt)
            obs=np.mean([np.abs(x@y) for x,y in zip(xs,ys)],axis=0); null=[]
            for _ in range(199):
                cc=[]
                for x,y,dt in zip(xs,ys,dts):
                    k=int(60/dt); shift=rng.integers(k,len(y)-k+1);cc.append(np.abs(x@np.roll(y,shift)))
                null.append(np.max(np.mean(cc,axis=0)))
            rows.append([sc,gas,int(LAGS[obs.argmax()]),obs.max(),np.quantile(null,.95),(1+np.sum(np.asarray(null)>=obs.max()))/200])
    pd.DataFrame(rows,columns=['scenario','gas','peak_s','mean_abs_r','null_p95','shift_p']).to_csv(out/'lag_negative_control.csv',index=False,float_format='%.8g')


def audit(out):
    out=Path(out); f=pd.read_csv(out/'fold_metrics.csv'); p=pd.read_csv(out/'predictions.csv.gz'); l=pd.read_csv(out/'lags.csv')
    keys=['scenario','gas','heat','model']; assert len(f)==760 and len(l)==40 and not f.duplicated(keys).any()
    assert not p.duplicated(keys+['detector_s']).any() and (p.prediction>=0).all()
    for _,z in p.groupby(['scenario','gas','heat']):
        y=z.pivot(index=['publish_s','detector_s'],columns='model',values='observed'); assert y.shape[1]==19 and y.notna().all().all() and (y.max(1)==y.min(1)).all()
    for k,z in p.groupby(keys):
        y=z.observed.to_numpy(); v=z.prediction.to_numpy(); err=v-y; m=f.set_index(keys).loc[k]
        assert len(z)==m.n and np.allclose(z.detector_s-z.publish_s,m.lead_s)
        np.testing.assert_allclose([np.abs(err).mean(),np.sqrt(np.mean(err**2)),1-np.sum(err**2)/np.sum((y-y.mean())**2)],m[['mae','rmse','r2']].to_numpy(float),rtol=1e-5,atol=1e-6)
    rng=np.random.default_rng(20261127); rows=[]
    comparisons=[('process_nonlinear','fusion'),('fusion','state'),('fusion','alpha_beta'),('fusion','soft_moe2'),('fusion','individual_ensemble'),('process_nonlinear_ensemble','individual_ensemble'),('process_nonlinear_ensemble','fusion_peak_ensemble'),('zero_alignment','fusion')]
    for (sc,gas),a in f.groupby(['scenario','gas']):
        m=a.pivot(index='heat',columns='model',values='mae')
        for b,c in comparisons:
            d=(m[b]-m[c]).to_numpy(); samples=rng.choice(d,(10000,len(d)),replace=True).mean(1); lo,hi=np.quantile(samples,[.025,.975]);rows.append([sc,gas,b,c,d.mean(),(d>0).sum(),lo,hi])
    pd.DataFrame(rows,columns=['scenario','gas','baseline','candidate','mae_gain','heats_improved','descriptive_ci_low','descriptive_ci_high']).to_csv(out/'paired_effects.csv',index=False,float_format='%.8g')
    sel=f[f.model!='zero_alignment'].sort_values(['inner_mae','model']).groupby(['scenario','gas','heat']).first().reset_index()
    sel.to_csv(out/'nested_family_selection.csv',index=False,float_format='%.8g')
    f.pivot(index=['scenario','gas','model'],columns='heat',values='mae').reset_index().to_csv(out/'fold_mae_matrix.csv',index=False,float_format='%.6f')
    print(f'PASS: {len(f)} fold metrics and {len(p)} predictions; equal targets, nonnegative outputs, and independent MAE/RMSE/R2 recomputation')
    print(sel.groupby(['scenario','gas']).agg(mae=('mae','mean'),rmse=('rmse','mean'),r2=('r2','mean'),worst=('mae','max')))
    return f,p


def figures(out,p):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dest=out/'figures';dest.mkdir(exist_ok=True)
    for (sc,gas),a in p.groupby(['scenario','gas']):
        fig,axes=plt.subplots(5,2,figsize=(12,12),sharey=True);hs=sorted(a.heat.unique())
        for ax,h in zip(axes.ravel(),hs):
            z=a[a.heat==h];v=z[z.model=='fusion_peak_ensemble'];w=z[z.model=='process_nonlinear_ensemble'];ax.plot(v.detector_s,v.observed,color='black',lw=.8,label='Observed');ax.plot(v.detector_s,v.prediction,lw=.9,label='Spectral ensemble');ax.plot(w.detector_s,w.prediction,lw=.8,label='Process ensemble');ax.set_title(h);ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=8);fig.supxlabel('Detector time (s)');fig.supylabel(f'{gas} flow (m3/s)');fig.suptitle(f'{sc} | Ten held-out heats');fig.tight_layout();fig.savefig(dest/f'{sc}_{gas}_predictions.png',dpi=150);plt.close(fig)


def main(raw,out):
    out=Path(out); f,p=audit(out); lag_control(raw,out);figures(out,p)
