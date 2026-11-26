"""Rebuild the optimized Q2 experiment and validate the generated output."""
import pandas as pd

from apmcm.common.paths import PROCESS_INPUT
from .checks import checks
from .delay import delay_study, local_profiles
from .figures import dose_response, figures
from .global_lag import full_profiles
from .reconstruct import ROOT, TAB, reconstruct
from .validate import validate


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    (ROOT / "figures").mkdir(exist_ok=True)
    curves = reconstruct(pd.read_csv(PROCESS_INPUT))
    delay_study(curves)
    local_profiles(curves)
    full_profiles(curves)
    checks(curves)
    figures()
    dose_response()
    validate()


if __name__ == "__main__":
    main()
