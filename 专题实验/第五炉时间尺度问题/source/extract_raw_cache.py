#!/usr/bin/env python3
"""Extract ten XLSB workbooks into compact NPZ caches without modifying source files."""
from pathlib import Path
import argparse, zipfile, re, subprocess, sys, json, hashlib, pandas as pd

p=argparse.ArgumentParser(); p.add_argument('--input',required=True,help='Data.zip or directory with ten XLSB files')
p.add_argument('--cache-dir',default='cache'); p.add_argument('--raw-dir',default='raw_data'); a=p.parse_args()
cache=Path(a.cache_dir).resolve(); cache.mkdir(parents=True,exist_ok=True)
src=Path(a.input).resolve()
if src.is_file() and src.suffix.lower()=='.zip':
 raw=Path(a.raw_dir).resolve(); raw.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(src) as z:
  members=[m for m in z.namelist() if m.lower().endswith('.xlsb')]
  if len(members)!=10: raise ValueError(f'Expected 10 XLSB workbooks; found {len(members)}')
  for m in members:
   dest=raw/Path(m).name
   if not dest.exists():
    with z.open(m) as r,dest.open('wb') as w:
     while True:
      b=r.read(1<<20)
      if not b: break
      w.write(b)
 files=list(raw.glob('*.xlsb'))
elif src.is_dir(): files=list(src.glob('*.xlsb'))
else: raise FileNotFoundError(src)
if len(files)!=10: raise ValueError(f'Expected 10 XLSB files; found {len(files)}')
files=sorted(files,key=lambda x:int(re.search(r'Attachment\s+(\d+)',x.name).group(1)))
worker=Path(__file__).with_name('extract_worker.py')
records=[]
for f in files:
 r=subprocess.run([sys.executable,str(worker),str(f),str(cache)],capture_output=True,text=True,check=True)
 rec=json.loads(r.stdout.strip()); print(rec); records.append(rec)
pd.DataFrame(records).to_csv(cache/'raw_manifest.csv',index=False)
if src.is_file():
 (cache/'DATA_SHA256.txt').write_text(f'Data.zip SHA-256: {hashlib.sha256(src.read_bytes()).hexdigest()}\n',encoding='utf-8')
print(f'Created {len(files)} NPZ caches in {cache}')
