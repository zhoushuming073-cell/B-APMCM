"""Render one Q1 surface in an isolated process."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

from .figures import _render_heat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--variant", choices=["hard_clip_lower_bound", "template_sensitivity"], required=True)
    parser.add_argument("--time-stride", type=int, required=True)
    parser.add_argument("--wavelength-stride", type=int, required=True)
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    args = parser.parse_args()

    fig, record, _ = _render_heat(
        args.data, args.variant, args.time_stride, args.wavelength_stride
    )
    args.image.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.image, dpi=120, bbox_inches="tight")
    plt.close(fig)
    args.record.write_text(json.dumps(record, ensure_ascii=False), encoding="utf-8")
    print(f"Rendered {record['heat']} ({args.variant})", flush=True)


if __name__ == "__main__":
    main()
