import argparse
from pathlib import Path

from apmcm.q2.specialized import run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reproduce this Q2 protocol")
    parser.add_argument("--output", type=Path)
    run("timebase_scenarios", parser.parse_args().output)
