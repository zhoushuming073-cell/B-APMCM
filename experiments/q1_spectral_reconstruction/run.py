import argparse
from pathlib import Path

from apmcm.common.paths import ROOT, config
from apmcm.q1.figures import draw
from apmcm.q1.reconstruct import reconstruct


def run_diagnostics():
    """Run the independent peak-masking and saturation validation studies."""
    from threadpoolctl import threadpool_limits
    from apmcm.q1.peak_study import study
    from apmcm.q1.peak_fit_metrics import main as peak_fit_main
    from apmcm.q1.saturation_study import main as saturation_main

    with threadpool_limits(limits=1):
        study()
        peak_fit_main()
        saturation_main()


if __name__ == "__main__":
    settings = config("q1.yaml")
    parser = argparse.ArgumentParser(
        description="Rebuild Q1 anomaly/saturation scenarios and native-grid spectral surfaces"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--heat", action="append", default=[])
    parser.add_argument("--time-stride", type=int, default=settings["display_time_stride"])
    parser.add_argument(
        "--wavelength-stride", type=int, default=settings["display_wavelength_stride"]
    )
    parser.add_argument(
        "--template-sensitivity", action="store_true",
        help="Also render the unconstrained temporal-template comparison PDF",
    )
    parser.add_argument(
        "--diagnostics", action="store_true",
        help="Also run peak masking, peak-fit, and saturation diagnostics",
    )
    args = parser.parse_args()

    output = args.output or ROOT / settings["rerun_output"]
    reconstruct(output, args.heat)
    draw(
        output / "scenarios",
        output / "figures",
        time_stride=args.time_stride,
        wavelength_stride=args.wavelength_stride,
        include_template_sensitivity=(
            settings["include_template_sensitivity"] or args.template_sensitivity
        ),
        heats=args.heat or None,
    )
    if args.diagnostics:
        run_diagnostics()
    print(f"Q1 finalization outputs written to: {output}")
