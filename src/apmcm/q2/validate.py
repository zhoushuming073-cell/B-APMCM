"""Check input identities, conservation, evaluation splits and main tables."""
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from .reconstruct import ROOT,TAB,build,CHANNELS
from apmcm.common.paths import PROCESS_INPUT


def validate():
    data=pd.read_csv(PROCESS_INPUT)
    assert len(data)==13589 and data.heat.nunique()==10
    assert data[['row','O','CO2','CO']].notna().all().all()
    assert (data.groupby('heat').O.diff().dropna()>=0).all()
    choices=pd.read_csv(TAB/'reconstruction_choice.csv').set_index(['scenario','heat','channel'])
    out=[]
    curves=pd.read_csv(TAB/'curves.csv.gz')
    for (s,h),g in curves.groupby(['scenario','heat']):
        t=g.t.to_numpy()
        dt=t[1]-t[0]
        assert len(g)==(data.heat==h).sum()
        assert g.t.is_monotonic_increasing
        assert g.q.min()>=-1.e-8 and np.diff(g.O).min()>=-1.e-8
        assert abs(g.O.iloc[0]-g.O_raw.iloc[0])<1.e-7 and abs(g.O.iloc[-1]-g.O_raw.iloc[-1])<1.e-7
        for v in ['E_CO2','E_CO']:
            assert g[v].iloc[0]==0 and np.diff(g[v]).min()>=-1.e-8
        for var in CHANNELS:
            method=choices.loc[(s,h,var),'method']
            f=build(t,g[var+'_raw'].to_numpy(),var,method,dt)
            assert np.max(abs(f(t)-g[var].to_numpy()))<1.e-5
            if var=='O':
                # Average-of-polygon derivative is piecewise linear between shifted knots.
                # Simpson on a fine grid also checks the PCHIP alternatives numerically.
                from scipy.integrate import simpson
                fine=np.linspace(t[0],t[-1],12*(len(t)-1)+1)
                err=abs(simpson(f.derivative()(fine),x=fine)-(f(t[-1])-f(t[0])))
                assert err<.01
                out.append(dict(scenario=s,heat=h,q_mass_error_m3=err))
    df=pd.read_csv(TAB/'delay_folds.csv')
    assert not df.duplicated(['scenario','window','heat','gas','model']).any()
    pred=pd.read_csv(TAB/'delay_predictions.csv.gz')
    for (s,h,gas),g in pred.groupby(['scenario','heat','gas']):
        piv=g.pivot(index='t',columns='model',values='raw')
        assert piv.shape[1]==5 and piv.notna().all().all() and (piv.max(axis=1)==piv.min(axis=1)).all()
        assert piv.index[0]==120
        for m,sub in g.groupby('model'):
            value=df[(df.scenario==s)&(df.heat==h)&(df.gas==gas)&(df.model==m)&(df.window==8)].mae.item()
            assert abs(abs(sub.raw-sub.estimated).mean()-value)<1.e-6
    assert len(list((ROOT/'figures').glob('heat_*.png')))==10
    pd.DataFrame(out).to_csv(TAB/'validation_mass.csv',index=False)
    print('PASS: 10 heats, 13,589 records; all 20 timebase/heat curves and 5-model paired held-out targets checked.')

if __name__=='__main__':
    validate()
