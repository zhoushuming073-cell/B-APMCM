import argparse
from pathlib import Path

from apmcm.common.paths import ROOT, config
from apmcm.q1.figures import draw
from apmcm.q1.reconstruct import reconstruct


if __name__ == "__main__":
    settings = config("q1.yaml")
    parser = argparse.ArgumentParser(description="Reconstruct and plot the original Q1 scenario model")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--heat", action="append", default=[])
    args = parser.parse_args()
    default = ROOT / settings["rerun_output"]
    if args.heat:
        default = default / "selected" / "_".join(sorted(set(args.heat)))
    output = args.output or default
    reconstruct(output, args.heat)
    draw(output / "scenarios", output / "figures", settings["display_stride"])
