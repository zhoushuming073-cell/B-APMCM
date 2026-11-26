import argparse
import os
from pathlib import Path

from apmcm.common.paths import ROOT, config


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reproduce the optimized Q2 experiment")
    parser.add_argument("--output", type=Path, default=ROOT / config("q2/optimized.yaml")["rerun_output"])
    os.environ["APMCM_Q2_OUTPUT"] = str(parser.parse_args().output.resolve())
    from apmcm.q2.__main__ import main
    main()
