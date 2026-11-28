"""Repository paths and the official input configuration."""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
CONFIG = ROOT / "configs"


def config(name):
    with (CONFIG / name).open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


official = config("data.yaml")
ARCHIVE = ROOT / official["archive"]
RAW = ROOT / official["workbooks"]
PROC = ROOT / official["spectral_cache"]
PROBLEM = ROOT / official["problem_statement"]
PROCESS_INPUT = ROOT / official["process_input"]
SOURCE_SHA256 = official["archive_sha256"]
DT = official["sample_interval_s"]
