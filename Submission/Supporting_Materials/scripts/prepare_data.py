"""Packaging-only helper: stage official inputs locally and rebuild extracted caches.
Never ship data/raw or data/interim/extracted in the supporting ZIP.
"""
import argparse
from pathlib import Path
from shutil import copyfile
from zipfile import ZipFile
import numpy as np
import pandas as pd

from apmcm.common.paths import ARCHIVE, RAW, PROC, DATA, PROCESS_INPUT, SOURCE_SHA256
from apmcm.data.verify import digest, verify
from apmcm.data.loader import run

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_zip", type=Path)
    args = parser.parse_args()
    src = args.data_zip.resolve()
    if digest(src) != SOURCE_SHA256:
        raise ValueError("Incorrect official Data.zip SHA-256")
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    if src != ARCHIVE.resolve():
        copyfile(src, ARCHIVE)
    with ZipFile(ARCHIVE) as z:
        entries = [n for n in z.namelist() if n.lower().endswith(".xlsb")]
        names = [Path(n).name for n in entries]
        if len(entries) != 10 or len(set(names)) != 10:
            raise ValueError("Expected 10 distinct workbook names")
        for entry, name in zip(entries, names):
            (RAW / name).write_bytes(z.read(entry))
    verify()
    run()
    manifest = pd.read_csv(DATA / "metadata/manifest.csv").sort_values("attachment")
    frames = []
    for item in manifest.itertuples(index=False):
        with np.load(PROC / f"raw_{item.heat}.npz", allow_pickle=False) as a:
            frames.append(pd.DataFrame({
                "row": a["sequence"].astype(int),
                "O": a["oxygen"],
                "CO2": a["gas"][:, 0],
                "CO": a["gas"][:, 1],
                "attachment": int(item.attachment),
                "heat": str(item.heat),
            }))
    process = pd.concat(frames, ignore_index=True)
    if len(process) != 13589 or process["heat"].nunique() != 10:
        raise ValueError("Unexpected process table dimensions")
    PROCESS_INPUT.parent.mkdir(parents=True, exist_ok=True)
    process.to_csv(PROCESS_INPUT, index=False)
    print("Prepared 10 heats / 13,589 rows. Keep original inputs and raw caches out of the ZIP.")

if __name__ == "__main__":
    main()
