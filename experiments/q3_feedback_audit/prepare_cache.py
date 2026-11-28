"""Build the compact causal Q3 cache used by run_feedback_sweep.py."""
from pathlib import Path
import argparse
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from apmcm.q3.data import load_archive


def build(archive: Path, output: Path):
    raw = load_archive(archive)
    output.mkdir(parents=True, exist_ok=True)
    edges = np.linspace(345, 1015, 33)
    for heat, a in raw.items():
        x = a["x"].astype(float)
        wl = a["wl"].astype(float)
        bands = np.stack([
            np.log1p(np.mean(x[:, (wl >= lo) & (wl < hi)], axis=1))
            for lo, hi in zip(edges[:-1], edges[1:])
        ], axis=1)
        contrast = np.stack([
            np.log1p(np.quantile(x[:, (wl >= lo) & (wl < hi)], .95, axis=1))
            - np.log1p(np.median(x[:, (wl >= lo) & (wl < hi)], axis=1))
            for lo, hi in zip(edges[:-1], edges[1:])
        ], axis=1)
        np.savez_compressed(
            output / f"{heat}.npz",
            O=a["oxygen"], bands=bands, contrast=contrast,
            gas=a["gas"], wl=wl,
        )
        print(heat, len(x), bands.shape, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--archive", type=Path, default=Path("data/raw/official/Data.zip"))
    p.add_argument("--out", type=Path, default=Path("data/interim/q3_feedback_cache"))
    args = p.parse_args()
    build(args.archive, args.out)
