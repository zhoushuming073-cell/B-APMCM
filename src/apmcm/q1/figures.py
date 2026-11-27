"""Render Q1 spectral surfaces at full native wavelength resolution."""
from __future__ import annotations

from pathlib import Path
import json
import os
import subprocess
import sys

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from apmcm.common.paths import DT
from .reconstruct import SATURATION_CEILING


def display_indices(length: int, stride: int) -> np.ndarray:
    """Select display nodes while retaining the final coordinate."""
    if length < 2:
        raise ValueError("A display axis must contain at least two nodes")
    if stride < 1:
        raise ValueError("Display stride must be a positive integer")
    return np.unique(np.r_[np.arange(0, length, stride), length - 1]).astype(int)


def _load_surface(data, variant: str) -> np.ndarray:
    if variant == "hard_clip_lower_bound":
        if "scenario_hard_clip_lower_bound" in data.files:
            return data["scenario_hard_clip_lower_bound"]
        surface = data["scenario"].copy()
        surface[data["saturated"]] = np.maximum(
            SATURATION_CEILING, surface[data["saturated"]]
        )
        return surface
    if variant == "template_sensitivity":
        return data["scenario_template"] if "scenario_template" in data.files else data["scenario"]
    raise ValueError(f"Unknown Q1 surface variant: {variant}")


def _render_heat(data_path: Path, variant: str, time_stride: int, wavelength_stride: int):
    heat = data_path.stem.removesuffix("_scenario")
    with np.load(data_path) as data:
        sequence = data["seq"].astype(float)
        wavelength = data["wl"].astype(float)
        surface = _load_surface(data, variant).astype(np.float32)
        saturated = data["saturated"].astype(bool)

    time = (sequence - sequence[0]) * DT
    time_rows = display_indices(len(time), time_stride)
    wavelength_cols = display_indices(len(wavelength), wavelength_stride)
    tt, ww = np.meshgrid(time[time_rows], wavelength[wavelength_cols], indexing="ij")
    zz = surface[np.ix_(time_rows, wavelength_cols)]

    fig = plt.figure(figsize=(12, 7.2))
    ax = fig.add_subplot(projection="3d")
    plotted = ax.plot_surface(
        tt, ww, zz, cmap="viridis", linewidth=0, antialiased=False,
        rstride=1, cstride=1,
    )
    sat_i, sat_j = np.where(saturated)
    if len(sat_i):
        ax.scatter(
            time[sat_i], wavelength[sat_j],
            np.full(len(sat_i), SATURATION_CEILING),
            color="crimson", marker=".", s=3, alpha=0.7,
            depthshade=False, rasterized=True,
            label=f"Recorded at ceiling (65,535): {len(sat_i):,} points",
        )
        ax.legend(loc="upper left", fontsize=8)

    note = ("conditional hard-clipping lower-bound scenario"
            if variant == "hard_clip_lower_bound"
            else "unconstrained temporal-template sensitivity scenario")
    suffix = "" if variant == "hard_clip_lower_bound" else "_template_sensitivity"
    ax.set(
        xlabel="Time (s)", ylabel="Wavelength (nm)",
        zlabel="Intensity (W/cm²)", title=f"{heat} | {note}",
    )
    ax.view_init(elev=27, azim=-135)
    fig.colorbar(plotted, ax=ax, shrink=0.58, pad=0.08, label="Intensity (W/cm²)")
    fig.tight_layout()
    record = {
        "heat": heat, "variant": variant,
        "time_rows_total": len(time), "time_rows_displayed": len(time_rows),
        "time_stride": time_stride,
        "wavelength_channels_total": len(wavelength),
        "wavelength_channels_displayed": len(wavelength_cols),
        "wavelength_stride": wavelength_stride,
        "saturation_points_marked": int(saturated.sum()),
        "time_display_interval_s": float(time[time_rows[1]] - time[time_rows[0]]),
    }
    return fig, record, suffix


def _write_image_pdf(image_paths: list[Path], pdf_path: Path, dpi: int = 120) -> None:
    """Create a compact PDF from rendered PNG pages, avoiding huge vector meshes."""
    document = fitz.open()
    try:
        for image_path in image_paths:
            pixmap = fitz.Pixmap(str(image_path))
            page = document.new_page(
                width=pixmap.width * 72.0 / dpi,
                height=pixmap.height * 72.0 / dpi,
            )
            page.insert_image(page.rect, filename=str(image_path))
        document.save(str(pdf_path), garbage=4, deflate=True)
    finally:
        document.close()


def draw(scenarios: Path, destination: Path, time_stride: int = 10,
         wavelength_stride: int = 1, include_template_sensitivity: bool = False,
         heats=None) -> pd.DataFrame:
    """Render each heat in a fresh process and assemble one PDF per variant.

    Time rows are decimated for display only. Wavelength channels default to
    stride 1, preserving all 2,047 native spectral channels.
    """
    scenarios, destination = Path(scenarios), Path(destination)
    all_files = sorted(scenarios.glob("*_scenario.npz"))
    if not all_files:
        raise FileNotFoundError(f"No reconstructed scenarios in {scenarios}")
    available = {p.stem.removesuffix("_scenario"): p for p in all_files}
    selected = sorted(set(heats)) if heats else sorted(available)
    missing = set(selected) - set(available)
    if missing:
        raise ValueError(f"Requested plot heats are missing: {sorted(missing)}")
    files = [available[h] for h in selected]
    full_selection = set(selected) == set(available)
    destination.mkdir(parents=True, exist_ok=True)

    variants = ["hard_clip_lower_bound"]
    if include_template_sensitivity:
        variants.append("template_sensitivity")
    records = []
    for variant in variants:
        suffix = "" if variant == "hard_clip_lower_bound" else "_template_sensitivity"
        image_paths = []
        for file in files:
            heat = file.stem.removesuffix("_scenario")
            image_path = destination / f"{heat}_spectral_surface{suffix}.png"
            record_path = destination / f".{heat}_{variant}.json"
            command = [
                sys.executable, "-m", "apmcm.q1.render_one",
                "--data", str(file), "--variant", variant,
                "--time-stride", str(time_stride),
                "--wavelength-stride", str(wavelength_stride),
                "--image", str(image_path), "--record", str(record_path),
            ]
            environment = os.environ.copy()
            environment["MPLCONFIGDIR"] = str(destination / f".mplconfig_{heat}_{variant}")
            subprocess.run(command, check=True, env=environment)
            records.append(json.loads(record_path.read_text(encoding="utf-8")))
            record_path.unlink(missing_ok=True)
            image_paths.append(image_path)

        pdf_stem = "all_heats" if full_selection else "selected_heats"
        _write_image_pdf(
            image_paths, destination / f"{pdf_stem}_spectral_surface{suffix}.pdf"
        )

    manifest = pd.DataFrame(records)
    manifest_name = "q1_plot_manifest.csv" if full_selection else "q1_plot_manifest_selected.csv"
    manifest.to_csv(destination / manifest_name, index=False)
    return manifest
