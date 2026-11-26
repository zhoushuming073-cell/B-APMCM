"""Positive controls, sensitivity, and uncertainty summarized at the heat level."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from reconstruct import TAB,GASES,flags,CHANNELS
from delay import series,design,CONFIG


def checks(curves):
    rng=np.random.default_rng(20261126)
    control=[]
    for heat,g in curves[curves.scenario=='H0'].groupby('heat'):
        for var in CHANNELS:
            y=g[var+'_raw'].to_numpy()
            t=g.t.to_numpy()
            baseline=flags(t,y,var,.5)
            positions=np.arange(15,len(y)-15,30)
            for amp in [1.,4.]:
                if var=='O': amplitude=15*amp
                else: amplitude=amp
                yc=y.copy()
                yc[positions]+=amplitude*rng.choice([-1,1],len(positions))
                flagged=flags(t,yc,var,.5)
                control.append(dict(heat=heat,channel=var,amplitude=amplitude,n=len(positions),
                                    recovered=flagged[positions].sum(),natural_candidates=baseline.sum(),
                                    additional_elsewhere=np.sum((flagged&~baseline)[np.setdiff1d(np.arange(len(y)),positions)])))
    pd.DataFrame(control).to_csv(TAB/'spike_positive_control.csv',index=False)
    # Positive control for the inversion and grid search using the observed oxygen inputs.
    # Known 24 s reference, two laws. Response has an exactly known linear q mapping here.
    trials=[]
    for truth in [0.,np.log(2.)]:
        costs=np.zeros((len(CONFIG)-1,40))
        for heat,g in curves[curves.scenario=='H0'].groupby('heat'):
            s=series(g,8)
            x,_=design(s,('shift',24.,truth))
            noise=rng.normal(0,.8,(len(x),40))
            target=4*x[:,-2,None]+noise
            for i,c in enumerate(CONFIG[1:]):
                z,_=design(s,c)
                costs[i]+=np.mean((4*z[:,-2,None]-target)**2,axis=0)
        for rep,i in enumerate(np.argmin(costs,axis=0)):
            _,d,b=CONFIG[i+1]
            trials.append(dict(truth_d=24.,truth_beta=truth,rep=rep+1,found_d=d,found_beta=b))
    pd.DataFrame(trials).to_csv(TAB/'delay_positive_control.csv',index=False)
    local=pd.read_csv(TAB/'local_windows.csv')
    slopes=[]
    for (scen,gas),group in local.groupby(['scenario','gas']):
        terms=[]
        for heat,g in group.groupby('heat'):
            time=g.t.to_numpy()
            X=np.column_stack([np.ones(len(g)),time-time.mean()])
            z=g[['q','lag']].to_numpy()
            residual=z-X@np.linalg.lstsq(X,z,rcond=None)[0]
            x,y=residual.T
            terms.append([x@x,x@y])
        terms=np.array(terms)
        idx=rng.integers(0,len(terms),(5000,len(terms)))
        b=terms[idx].sum(axis=1)
        lo,hi=np.quantile(b[:,1]/b[:,0],[.025,.975])
        slopes.append(dict(scenario=scen,gas=gas,n_heats=len(terms),slope=terms[:,1].sum()/terms[:,0].sum(),
                           lo=lo,hi=hi,bootstrap=5000))
    pd.DataFrame(slopes).to_csv(TAB/'flow_slope.csv',index=False)
    folds=pd.read_csv(TAB/'delay_folds.csv')
    paired=[]
    for (scen,w,gas),g in folds.groupby(['scenario','window','gas']):
        pivot=g.pivot(index='heat',columns='model',values='mae')
        diff=(pivot.decreasing-pivot.fixed).to_numpy()
        means=diff[rng.integers(0,len(diff),(5000,len(diff)))].mean(axis=1)
        paired.append(dict(scenario=scen,window=w,gas=gas,delta_mae=diff.mean(),lo=np.quantile(means,.025),
                           hi=np.quantile(means,.975),improved=(diff < -1.e-8).sum(),tied=(abs(diff)<=1.e-8).sum(),
                           worse=(diff>1.e-8).sum()))
    pd.DataFrame(paired).to_csv(TAB/'paired_delay.csv',index=False)
    audit=[]
    cv=pd.read_csv(TAB/'reconstruction_cv.csv')
    choices=pd.read_csv(TAB/'reconstruction_choice.csv')
    for label,df in [('selected',cv.merge(choices)),('linear baseline',cv[cv.method=='linear'])]:
        for (scen,heat,var),g in df.groupby(['scenario','heat','channel']):
            sse=(g.n*g.rmse**2).sum()
            n=g.n.sum()
            r2=1-sse/(g.sum_y2.sum()-g.sum_y.sum()**2/n)
            audit.append(dict(scenario=scen,heat=heat,channel=var,model=label,n=n,
                              mae=(g.n*g.mae).sum()/n,rmse=np.sqrt(sse/n),r2=r2))
    pd.DataFrame(audit).to_csv(TAB/'reconstruction_holdout.csv',index=False)
    totals=[]
    for (scen,heat),g in curves.groupby(['scenario','heat']):
        row=dict(scenario=scen,heat=heat,duration=g.t.iloc[-1],net_O=g.O.iloc[-1]-g.O.iloc[0])
        for gas in GASES:
            row['E_'+gas]=g['E_'+gas].iloc[-1]
            row['raw_trap_'+gas]=np.trapezoid(g[gas+'_raw'],g.t)
        totals.append(row)
    pd.DataFrame(totals).to_csv(TAB/'emissions.csv',index=False)
