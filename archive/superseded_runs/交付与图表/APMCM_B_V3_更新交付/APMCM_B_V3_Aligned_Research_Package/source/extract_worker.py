#!/usr/bin/env python3
"""One-workbook XLSB reader; run as a separate process for robust resource isolation."""
from pathlib import Path
import sys, re, hashlib, json
import numpy as np
import pandas as pd

src=Path(sys.argv[1]); cache_dir=Path(sys.argv[2]); cache_dir.mkdir(parents=True,exist_ok=True)
m=re.search(r'Heat No\.\s*(.+)\.xlsb$',src.name)
if not m: raise ValueError(f'Unrecognized workbook name: {src.name}')
heat=m.group(1)
df=pd.read_excel(src,engine='pyxlsb')
arr=df.to_numpy(dtype=np.float64)
if arr.ndim!=2 or arr.shape[1]!=2051: raise ValueError(f'{src.name}: expected 2051 columns, got {arr.shape}')
seq=arr[:,0]; X=arr[:,1:2048]; proc=arr[:,2048:2051]
wl=np.asarray(df.columns[1:2048],dtype=np.float64)
if len(wl)!=2047: raise ValueError(f'{src.name}: expected 2047 wavelengths, got {len(wl)}')
out=cache_dir/f'heats_{heat}.npz'
np.savez_compressed(out,seq=seq,X=X,proc=proc,wl=wl)
print(json.dumps({'heat':heat,'shape':list(arr.shape),'cache':str(out),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
       'sequence_contiguous':bool(np.array_equal(seq,np.arange(1,len(seq)+1))),
       'missing_cells':int(np.isnan(arr).sum()),'duplicate_sequence':int(pd.Series(seq).duplicated().sum()),
       'duplicate_rows_excluding_sequence':int(pd.DataFrame(arr[:,1:]).duplicated().sum())}))
