"""Generate Q1 surfaces and the optimized Q2 figures from published results."""
import argparse
from pathlib import Path

from apmcm.common.paths import RESULTS, ROOT, config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--question", choices=("q1", "q2", "all"), default="all")
    parser.add_argument("--output", type=Path, default=RESULTS / "rebuilt_figures")
    args = parser.parse_args()
    if args.question in ("q1", "all"):
        from apmcm.q1.figures import draw
        settings = config("q1.yaml")
        draw(ROOT / settings["published_scenarios"], args.output / "q1",\n             time_stride=settings["display_time_stride"],\n             wavelength_stride=settings["display_wavelength_stride"])
    if args.question in ("q2", "all"):
        from apmcm.q2 import figures
        figures.TAB = ROOT / config("q2/optimized.yaml")["published_tables"]
        figures.OUT = args.output / "q2"
        figures.OUT.mkdir(parents=True, exist_ok=True)
        figures.figures()


if __name__ == "__main__":
    main()
