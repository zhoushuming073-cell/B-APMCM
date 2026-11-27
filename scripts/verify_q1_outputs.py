"""Independently validate saved Q1 scenarios and rebuild aggregate audit tables."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from apmcm.common.paths import DT, PROC, ROOT, SOURCE_SHA256, config

CEILING = 65535.0


def verify(output: Path) -> pd.DataFrame:
    scenarios_dir = output / "scenarios"
    tables_dir = output / "tables"
    files = sorted(scenarios_dir.glob("*_scenario.npz"))
    expected = {p.stem[4:] for p in PROC.glob("raw_*.npz")}
    actual = {p.stem.removesuffix("_scenario") for p in files}
    if len(expected) != 10 or actual != expected:
        raise ValueError(
            f"Expected all ten official heats; missing={sorted(expected-actual)}, "
            f"unexpected={sorted(actual-expected)}"
        )

    summary, saturation_rows = [], []
    for path in files:
        heat = path.stem.removesuffix("_scenario")
        with np.load(path) as data:
            required = {
                "raw", "scenario", "scenario_anomaly_only", "scenario_template",
                "scenario_hard_clip_lower_bound", "wl", "seq", "score", "screen",
                "approved", "saturated", "saturation_repaired",
                "saturation_template_values", "saturation_template_below_ceiling",
                "saturation_template_above_ceiling", "saturation_unresolved",
                "saturation_ceiling", "source_sha256", "detector_threshold",
            }
            missing = required - set(data.files)
            if missing:
                raise ValueError(f"{heat}: missing arrays {sorted(missing)}")
            raw = data["raw"]
            anomaly = data["scenario_anomaly_only"]
            template = data["scenario_template"]
            legacy = data["scenario"]
            lower = data["scenario_hard_clip_lower_bound"]
            wavelength, sequence, score = data["wl"], data["seq"], data["score"]
            screen, approved = data["screen"].astype(bool), data["approved"].astype(bool)
            saturated = data["saturated"].astype(bool)
            repaired = data["saturation_repaired"].astype(bool)
            template_values = data["saturation_template_values"]
            below = data["saturation_template_below_ceiling"].astype(bool)
            above = data["saturation_template_above_ceiling"].astype(bool)
            unresolved = data["saturation_unresolved"].astype(bool)
            ceiling = float(data["saturation_ceiling"])
            source_hash = str(data["source_sha256"])
            threshold = float(data["detector_threshold"])

        if source_hash != SOURCE_SHA256 or ceiling != CEILING:
            raise ValueError(f"{heat}: source hash or saturation ceiling mismatch")
        if raw.ndim != 2 or raw.shape[1] != 2047 or raw.shape[0] != len(sequence):
            raise ValueError(f"{heat}: unexpected raw matrix shape")
        if len(wavelength) != 2047 or not np.all(np.diff(wavelength) > 0):
            raise ValueError(f"{heat}: invalid native wavelength grid")
        if not np.all(np.diff(sequence) == 1):
            raise ValueError(f"{heat}: sequence index is discontinuous")
        arrays = (anomaly, template, legacy, lower, score, screen, approved,
                  saturated, repaired, template_values, below, above, unresolved)
        if any(a.shape != raw.shape for a in arrays):
            raise ValueError(f"{heat}: output array shape mismatch")
        if not np.isfinite(score).all() or not np.isfinite(template).all():
            raise ValueError(f"{heat}: non-finite score or template scenario")
        if not np.array_equal(saturated, raw >= ceiling):
            raise ValueError(f"{heat}: saturation mask differs from raw observations")
        if (approved & ~screen).any() or (screen & saturated).any():
            raise ValueError(f"{heat}: screening/approval masks are inconsistent")
        if (repaired & ~saturated).any():
            raise ValueError(f"{heat}: saturation estimate outside saturation mask")
        if not np.array_equal(legacy, template):
            raise ValueError(f"{heat}: legacy scenario alias differs from template scenario")
        if not np.array_equal(anomaly[~approved], raw[~approved]):
            raise ValueError(f"{heat}: anomaly scenario changed an unapproved point")
        changed = approved | repaired
        if not np.array_equal(template[~changed], raw[~changed]):
            raise ValueError(f"{heat}: template scenario changed a point outside its masks")
        if not np.isfinite(template_values[repaired]).all() or not np.isnan(template_values[~repaired]).all():
            raise ValueError(f"{heat}: saturation template values do not match their mask")
        if not np.array_equal(below, saturated & repaired & (template_values < ceiling)):
            raise ValueError(f"{heat}: below-ceiling mask is inconsistent")
        if not np.array_equal(above, saturated & repaired & ~below):
            raise ValueError(f"{heat}: above-ceiling mask is inconsistent")
        if not np.array_equal(unresolved, saturated & ~repaired):
            raise ValueError(f"{heat}: unresolved mask is inconsistent")
        if (lower[saturated] < ceiling).any():
            raise ValueError(f"{heat}: hard-clipping lower bound falls below ceiling")
        if not np.array_equal(lower[~saturated], template[~saturated]):
            raise ValueError(f"{heat}: lower-bound scenario changed an unsaturated point")

        summary.append({
            "heat": heat, "n_time": raw.shape[0], "n_wavelength": raw.shape[1],
            "screen_candidates": int(screen.sum()), "screen_rate": float(screen.mean()),
            "approved_anomaly_repairs": int(approved.sum()), "approved_rate": float(approved.mean()),
            "saturation_points": int(saturated.sum()), "saturation_repaired": int(repaired.sum()),
            "saturation_template_below_ceiling": int(below.sum()),
            "saturation_template_above_ceiling": int(above.sum()),
            "saturation_unresolved": int(unresolved.sum()),
            "total_changed": int(approved.sum() + repaired.sum()),
            "threshold": threshold, "source_sha256": source_hash,
        })
        si, sj = np.where(saturated)
        if len(si):
            saturation_rows.append(pd.DataFrame({
                "heat": heat, "sequence": sequence[si],
                "time_s": (sequence[si] - sequence[0]) * DT,
                "wavelength_nm": wavelength[sj], "raw": raw[si, sj],
                "template_estimate": template_values[si, sj],
                "template_available": repaired[si, sj],
                "template_below_ceiling": below[si, sj],
                "hard_clip_lower_bound": lower[si, sj],
            }))

    result = pd.DataFrame(summary).sort_values("heat").reset_index(drop=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(tables_dir / "q1_full_reconstruction_summary.csv", index=False)
    audit = pd.concat(saturation_rows, ignore_index=True) if saturation_rows else pd.DataFrame()
    audit.to_csv(tables_dir / "q1_saturation_audit.csv", index=False)
    return result


def main():
    settings = config("q1.yaml")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / settings["rerun_output"])
    args = parser.parse_args()
    result = verify(args.output)
    columns = [
        "screen_candidates", "approved_anomaly_repairs", "saturation_points",
        "saturation_repaired", "saturation_template_below_ceiling", "saturation_unresolved",
    ]
    print(result.to_string(index=False))
    print("TOTAL", result[columns].sum().to_dict())
    print(f"Verified ten heat scenarios in {args.output}")


if __name__ == "__main__":
    main()
