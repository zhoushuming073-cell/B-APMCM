"""Draw ten Q1 scenario surfaces from the saved full-resolution grids."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from apmcm.common.paths import DT


def draw(scenarios, destination, stride=10):
    files = sorted(scenarios.glob("*_scenario.npz"))
    if not files:
        raise FileNotFoundError(f"No reconstructed scenarios in {scenarios}")
    destination.mkdir(parents=True, exist_ok=True)
    with PdfPages(destination / "all_heats_spectral_surface.pdf") as pdf:
        for file in files:
            heat = file.stem.removesuffix("_scenario")
            with np.load(file) as data:
                seq = data["seq"]
                t = (seq - seq[0]) * DT
                w, z = data["wl"], data["scenario"]
            rows = np.unique(np.r_[np.arange(0, len(t), stride), len(t)-1])
            cols = np.unique(np.r_[np.arange(0, len(w), stride), len(w)-1])
            tt, ww = np.meshgrid(t[rows], w[cols], indexing="ij")
            fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": "3d"})
            ax.plot_surface(tt, ww, z[np.ix_(rows, cols)], cmap="viridis", linewidth=0,
                            rstride=1, cstride=1)
            ax.set(xlabel="Time (s)", ylabel="Wavelength (nm)", zlabel="Recorded intensity",
                   title=f"{heat} | reconstructed flame spectrum")
            ax.view_init(elev=27, azim=-135)
            fig.tight_layout()
            fig.savefig(destination / f"{heat}_spectral_surface.png", dpi=140)
            pdf.savefig(fig)
            plt.close(fig)
