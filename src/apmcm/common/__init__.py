from pathlib import Path
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .paths import ROOT, DATA, RESULTS, SOURCE_SHA256, RAW, PROC, DT
CURRENT_DATA = DATA / "metadata"
TABLE = RESULTS / "q1/rerun/baseline/tables"
FIG = RESULTS / "q1/rerun/baseline/figures"
MODEL = RESULTS / "q1/rerun/baseline/spectra"
REPORT = ROOT / "reports"
GASES = ["CO2", "CO"]
SEED = 2601

def setup():
    for folder in [TABLE, FIG, MODEL]:
        folder.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(TABLE.parent / "run.log", encoding="utf8")])
    plt.rcParams.update({"font.size": 9, "figure.dpi": 110, "axes.spines.top": False})

def table(name, records):
    frame = pd.DataFrame(records)
    frame.to_csv(TABLE / (name+".csv"), index=False)
    return frame

def figure(fig, name):
    fig.savefig(FIG / (name+".png"), dpi=145, bbox_inches="tight")
    plt.close(fig)

def loads():
    paths = sorted(PROC.glob("raw_*.npz"))
    if len(paths) != 10:
        raise ValueError("Corrected dataset is incomplete. Run python scripts/verify_data.py --rebuild-cache first.")
    result = {}
    for p in paths:
        with np.load(p) as values:
            result[p.stem.split("_")[1]] = dict(values)
    return result
