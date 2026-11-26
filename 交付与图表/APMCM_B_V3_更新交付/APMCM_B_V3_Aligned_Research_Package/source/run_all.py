#!/usr/bin/env python3
"""One-command orchestration for the time-base sensitivity study."""
from pathlib import Path
import argparse, subprocess, sys
p=argparse.ArgumentParser(); p.add_argument('--input',required=True,help='Data.zip or extracted XLSB directory')
p.add_argument('--workdir',default='apmcm_timebase_results'); a=p.parse_args()
root=Path(__file__).resolve().parent; work=Path(a.workdir).resolve(); cache=work/'cache'; raw=work/'raw_data'; out=work/'results'
commands=[
 [sys.executable,str(root/'extract_raw_cache.py'),'--input',str(Path(a.input).resolve()),'--cache-dir',str(cache),'--raw-dir',str(raw)],
 [sys.executable,str(root/'analyze_audit.py'),'--cache-dir',str(cache),'--out',str(out)],
 [sys.executable,str(root/'fit_lag_models.py'),'--cache-dir',str(cache),'--out',str(out)],
 [sys.executable,str(root/'make_figures.py'),'--cache-dir',str(cache),'--out',str(out)],
]
for cmd in commands:
 print('RUN',' '.join(cmd),flush=True); subprocess.run(cmd,check=True)
print('Finished. Results:',out)
