#!/usr/bin/env python3
"""Launcher for the archived Q3 stage-aware research runner."""
from pathlib import Path
import gzip

source = Path(__file__).with_name("stage_aware_source.py.gz")
namespace = {"__name__": "__main__", "__file__": str(source)}
exec(compile(gzip.decompress(source.read_bytes()), str(source), "exec"), namespace)
