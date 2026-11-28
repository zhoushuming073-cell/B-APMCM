#!/usr/bin/env python3
"""One-command Q2.3 final validation from the official Data.zip input."""
from __future__ import annotations
import argparse, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from extract import extract
from build_curves import build
import q23_final_validation as analysis
from audit_archives import run as audit_archives
from write_report import write_report

def main():
    p=argparse.ArgumentParser(description='Rebuild and validate B-APMCM Q2.3 from the official input archive.')
    p.add_argument('--data-zip',type=Path,default=ROOT/'Data.zip')
    p.add_argument('--out',type=Path,default=ROOT/'q23_final_output')
    p.add_argument('--simulations',type=int,default=100)
    p.add_argument('--seed',type=int,default=20261128)
    a=p.parse_args()
    archive=a.data_zip.resolve();out=a.out.resolve()
    if not archive.exists():raise FileNotFoundError(f'Official input archive not found: {archive}')
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise RuntimeError(f'Output directory must be empty: {out}')
    inputs=out/'inputs';curves_dir=out/'curves';results=out/'results'
    print(f'[1/3] Extracting official archive: {archive}',flush=True)
    extract(archive,inputs)
    print('[2/3] Rebuilding frozen Q2.1 continuous curves',flush=True)
    flow_profile=ROOT/'inputs'/'frozen_q2_1'/'flow_profile.csv'
    curves,manifest,qa=build(inputs/'process.csv',curves_dir,archive,flow_profile)
    print(f"Rebuilt {manifest['rows']} curve rows; SHA256={manifest['curve_csv_sha256']}",flush=True)
    print('[3/3] Running Q2.3 final validation and 100-replicate semi-synthetic experiments',flush=True)
    sys.argv=['q23_final_validation.py','--process',str(inputs/'process.csv'),'--curves',str(curves_dir/'final_curves.csv'),
              '--out',str(results),'--archive',str(archive),'--simulations',str(a.simulations),'--seed',str(a.seed)]
    analysis.main()
    print('[4/4] Comparing archived references and writing final report',flush=True)
    audit_archives(results,ROOT/'inputs'/'archive_reference',ROOT/'inputs'/'frozen_q2_1',ROOT/'inputs'/'archive_reference'/'round3_input_curves.csv')
    qa_path=results/'qa_manifest.json'
    import json
    qadata=json.loads(qa_path.read_text());qadata['flow_profile_sha256']=manifest['flow_profile_sha256'];qadata['flow_profile_source']=manifest['flow_profile_source'];qa_path.write_text(json.dumps(qadata,indent=2),encoding='utf-8')
    write_report(results,results/'Q2_3_FINAL_VALIDATION_AND_FREEZE.md')
    print(f'FINAL OUTPUT: {out}',flush=True)

if __name__=='__main__':main()
