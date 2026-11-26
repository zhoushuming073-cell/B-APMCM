from pathlib import Path
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CURRENT_DATA = ROOT / "data/corrected_2026_11_26"
SOURCE_SHA256 = "cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727"
RAW = CURRENT_DATA / "raw"
PROC = CURRENT_DATA / "processed"
TABLE = ROOT / "outputs/tables"
FIG = ROOT / "outputs/figures"
MODEL = ROOT / "outputs/models"
REPORT = ROOT / "reports"
DT = 0.5
GASES = ["CO2", "CO"]
SEED = 2601

def setup():
    for folder in [RAW, PROC, TABLE, FIG, MODEL, REPORT]:
        folder.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(ROOT/"outputs/run.log", encoding="utf8")])
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
        raise ValueError("Corrected dataset is incomplete. Run src.data_loader.run() first.")
    result = {}
    for p in paths:
        with np.load(p) as values:
            result[p.stem.split("_")[1]] = dict(values)
    return result
