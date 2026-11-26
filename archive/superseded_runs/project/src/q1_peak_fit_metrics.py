"""Score the existing seven-channel Q1 peak masking protocol on corrected spectra."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from sklearn.utils.extmath import randomized_svd
from threadpoolctl import threadpool_limits

PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "data/corrected_2026_11_26/processed"
REFERENCE = PROJECT.parent / "APMCM_Full_Experiments/APMCM_Full_Experiments/q1_native_basis_long.csv"
OUT = PROJECT / "outputs/tables"


def peak_locations(x, heat):
    saturated = (x >= 65535).any(axis=1)
    rng = np.random.default_rng(20261127 + sum(map(ord, heat)))
    candidates = []
    for i in range(25, len(x) - 25, 25):
        if saturated[i]:
            continue
        peaks, _ = find_peaks(x[i], prominence=max(100, np.median(x[i]) * .08),
                              width=(1, 20), distance=4)
        peaks = [int(j) for j in peaks if 60 <= j < x.shape[1] - 60]
        rng.shuffle(peaks)
        candidates.extend((i, j) for j in peaks[:2])
    rng.shuffle(candidates)
    return candidates[:60]


def scores(points):
    y = points.actual.to_numpy()
    p = points.predicted.to_numpy()
    errors = y - p
    # Use the same observable shoulder baseline for the truth and prediction.
    relative_y = y - points.shoulder.to_numpy()
    denom = np.sum((relative_y - relative_y.mean()) ** 2)
    return dict(n_points=len(y), pearson_r=np.corrcoef(y, p)[0, 1],
                r2=1 - errors @ errors / np.sum((y - y.mean()) ** 2),
                r2_shoulder_subtracted=1 - errors @ errors / denom,
                mae=np.mean(np.abs(errors)))


def main():
    ref = pd.read_csv(REFERENCE)
    heats = list(ref.heat.unique())  # Keep the source experiment's stacking order.
    spectra = {}
    training = {}
    for heat in heats:
        with np.load(SOURCE / f"raw_{heat}.npz") as data:
            x = data["x"].astype(np.float32)
        spectra[heat] = x
        ids = np.linspace(0, len(x)-1, min(80, len(x))).astype(int)
        ids = ids[~(x[ids] >= 65535).any(axis=1)]
        sample = np.log1p(x[ids])
        training[heat] = sample - np.median(sample, axis=1, keepdims=True)

    records = []
    for heat, x in spectra.items():
        train = np.vstack([training[h] for h in heats if h != heat and h[:2] == heat[:2]])
        mu = train.mean(axis=0)
        _, _, basis = randomized_svd(train - mu, n_components=12, n_iter=5, random_state=13)
        for i, j in peak_locations(x, heat):
            hidden = np.arange(j - 3, j + 4)
            visible = np.ones(x.shape[1], dtype=bool)
            visible[hidden] = False
            row = x[i]
            shoulder = np.median(row[np.r_[j-11:j-3, j+4:j+12]])
            y = np.log1p(row)
            level = np.median(y[visible])
            coef = np.linalg.lstsq(basis[:, visible].T, (y-level)[visible]-mu[visible],
                                    rcond=None)[0]
            pred = np.expm1(np.maximum(0, mu + coef @ basis + level))[hidden]
            linear = np.interp(hidden, np.flatnonzero(visible), row[visible])
            for method, values in (("same_native_r12", pred), ("linear", linear)):
                for c, estimate in zip(hidden, values):
                    records.append((heat, i, j, c, method, row[c], estimate, shoulder))

    points = pd.DataFrame(records, columns=["heat", "frame", "peak_channel", "channel",
                                            "method", "actual", "predicted", "shoulder"])
    differences = (points.actual - points.predicted).abs()
    event = points.assign(MAE=differences).groupby(
        ["heat", "method", "frame", "peak_channel"], sort=False).MAE.mean()
    for (heat, method), group in event.groupby(level=["heat", "method"]):
        reference = ref[(ref.heat == heat) & (ref.width == 7) & (ref.method == method)]
        # Check individual events in their original order before publishing any scores.
        np.testing.assert_allclose(group.to_numpy(), reference.MAE, rtol=3e-3, atol=0.2)
    result = [dict(heat=heat, method=method, **scores(group))
              for (heat, method), group in points.groupby(["heat", "method"])]
    result.extend(dict(heat="pooled", method=method, **scores(group))
                  for method, group in points.groupby("method"))
    pd.DataFrame(result).to_csv(OUT / "q1_peak_mask_fit_metrics.csv", index=False)
    points.to_csv(OUT / "q1_peak_mask_fit_points.csv", index=False)


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
