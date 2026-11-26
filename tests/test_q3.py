"""Regression checks for Q3 causal feature extraction and held-out labels."""
import numpy as np
from apmcm.q3.model import features, prepare, lag, predict_candidates


def synthetic():
    rng=np.random.default_rng(123)
    wl=np.linspace(339,1028,2047)
    raw={}
    for i in range(5):
        t=np.arange(320)*.5
        x=(200+100*np.sin(t[:,None]/15+i)+rng.uniform(0,100,(len(t),len(wl)))).astype('f4')
        raw[str(i)]=dict(wl=wl,x=x,oxygen=4*t,gas=np.c_[8+np.sin(t/20+i),20+4*np.cos(t/25+i)])
    return raw


def test_no_future_features():
    a=synthetic()['0']; j=160
    x=a['x'].copy();o=a['oxygen'].copy()
    x[j+1:]*=2;o[j+1:]+=1234
    f=features(a['wl'],a['x'],a['oxygen'],.5); g=features(a['wl'],x,o,.5)
    for i in (0,1,2,4):
        np.testing.assert_array_equal(f[i][:j+1],g[i][:j+1])
    np.testing.assert_array_equal(f[3],g[3])


def test_heldout_gas_not_used():
    d=prepare(synthetic(),'H0');test='4';train=['0','1','2','3'];lead=lag(d,train,0)
    a,_=predict_candidates(d,train,test,0,lead)
    d[test]['gas']+=10000;d[test]['profile'][:]=1
    assert lag(d,train,0)==lead
    b,_=predict_candidates(d,train,test,0,lead)
    for key in a:
        np.testing.assert_array_equal(a[key],b[key])
