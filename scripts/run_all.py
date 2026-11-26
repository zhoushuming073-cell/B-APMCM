"""Verify the data and run selected experiments in separate rerun directories."""
import argparse
import os
import subprocess
import sys

from apmcm.common.paths import ROOT

EXPERIMENTS = ('q1_spectral_reconstruction', 'q2_aligned_comparison', 'q2_delay_identification', 'q2_timebase_audit', 'q2_optimized')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiments", nargs="*", metavar="EXPERIMENT")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    selected = args.experiments or EXPERIMENTS
    if unknown := set(selected) - set(EXPERIMENTS):
        parser.error(f"Unknown experiments: {', '.join(sorted(unknown))}")
    subprocess.run([sys.executable, str(ROOT / "scripts/verify_data.py")], check=True)
    if args.verify_only:
        return
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    for name in selected:
        print(f"Running {name}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "experiments" / name / "run.py")],
                       cwd=ROOT, env=env, check=True)


if __name__ == "__main__":
    main()
