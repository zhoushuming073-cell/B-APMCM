"""Reproduce the B-only pilot from the official archive or exact bundled feature cache."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

from apmcm.common.paths import ROOT, config
from apmcm.q3.data import load_archive

HEATS = ('6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108')

def build_cache(archive, cache):
    raw = load_archive(archive)
    cache.mkdir(parents=True, exist_ok=True)
    edges = np.linspace(345, 1015, 33)
    for heat, a in raw.items():
        wl, x = a['wl'], a['x']
        bands, contrast = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            b = x[:, (wl >= lo) & (wl < hi)].astype(np.float64)
            bands.append(np.log1p(b.mean(axis=1)))
            contrast.append(np.log1p(np.quantile(b, .95, axis=1))-np.log1p(np.median(b, axis=1)))
        np.savez_compressed(cache/(heat+'.npz'), O=a['oxygen'].astype('f4'),
                            gas=a['gas'].astype('f4'), bands=np.stack(bands,axis=1).astype('f4'),
                            contrast=np.stack(contrast,axis=1).astype('f4'),
                            wl=wl, sat=(x>=65535).mean(axis=1))

def main():
    cfg=config('q3/b_holdout.yaml')
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,default=ROOT/cfg['archive'])
    p.add_argument('--cache',type=Path,help='Exact q3_eda directory from the original experiment ZIP (optional)')
    p.add_argument('--output',type=Path,default=ROOT/cfg['rerun_output'])
    p.add_argument('--versions',nargs='+',type=int,choices=(1,2,3),default=cfg['versions'])
    args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    cache=args.cache.resolve() if args.cache else out/'feature_cache'
    if args.cache is None:build_cache(args.archive,cache)
    if tuple(sorted(p.stem for p in cache.glob('6*.npz')))!=HEATS:
        raise ValueError('Expected exactly the official ten feature files')
    for version in args.versions:
        dest=out/f'v{version}';dest.mkdir(exist_ok=True)
        env=dict(os.environ,Q3_B_CACHE=str(cache),Q3_B_OUTPUT=str(dest),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,'-m',f'apmcm.q3.b_trial_v{version}'],cwd=ROOT,env=env,check=True)

if __name__=='__main__':main()
