#!/usr/bin/env python3
"""Repaired, reproducible second-round 40-test harness.

The archived GitHub entrypoint contained a one-character syntax error in the
AR(1) prewhitening slice (`x[1;]` -> `x[1:]`). The full 40-test family is
reproduced by the shared implementation in q23_final_validation.py.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from q23_final_validation import round2_40

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--curves',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    curves=pd.read_csv(a.curves)
    windows,summary=round2_40(curves,a.out)
    primary=summary[summary.offset_s==0]
    print(f'primary tests={len(primary)}; Holm q<=0.05={(primary.holm_p_40<=.05).sum()}')

if __name__=='__main__':main()
