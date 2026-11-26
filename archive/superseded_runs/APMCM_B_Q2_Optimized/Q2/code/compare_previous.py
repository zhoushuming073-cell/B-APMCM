"""Optional independent comparison against the input curves in the earlier experiment."""
import sys
from zipfile import ZipFile
import numpy as np
import pandas as pd
from reconstruct import ROOT,TAB

def compare(path):
    data=pd.read_csv(ROOT/'data/process.csv')
    rows=[]
    with ZipFile(path) as older:
        for heat,g in data.groupby('heat'):
            frame=pd.read_csv(older.open(f'inputs/curves_{heat}.csv'))
            assert len(frame)==len(g)
            for col,old in [('O','oxygen_cumulative'),('CO2','CO2_flow'),('CO','CO_flow')]:
                err=abs(g[col].to_numpy()-frame[old].to_numpy())
                rows.append(dict(heat=heat,channel=col,rows=len(g),max_abs=err.max()))
    tab=pd.DataFrame(rows)
    tab.to_csv(TAB/'raw_crosscheck.csv',index=False)
    assert tab.max_abs.max()==0

if __name__=='__main__':
    compare(sys.argv[1])
