"""Run Q3 from the canonical ZIP, then independently audit predictions."""
import argparse
from pathlib import Path
from threadpoolctl import threadpool_limits

from apmcm.common.paths import ROOT, config
from apmcm.q3.data import load_archive
from apmcm.q3.model import save
from apmcm.q3.diagnostics import main

if __name__ == '__main__':
    cfg=config('q3/spectral_gas_prediction.yaml')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ROOT/cfg['archive'])
    parser.add_argument('--output', type=Path, default=ROOT/cfg['rerun_output'])
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        raw=load_archive(args.archive)
        save(raw,args.output)
        main(raw,args.output)
