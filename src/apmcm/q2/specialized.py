"""Run the three archived Q2 protocols with canonical inputs and separate outputs."""
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

import pandas as pd

from apmcm.common.paths import ARCHIVE, ROOT, config


def run(protocol, output=None):
    settings = config(f"q2/{protocol}.yaml")
    source = ROOT / settings["frozen_source"]
    destination = output or ROOT / settings["rerun_output"]
    if protocol == "timebase_scenarios":
        entry = settings["scripts"][0]
        subprocess.run([sys.executable, "-X", "utf8", str(source / entry), "--input", str(ARCHIVE),
                        "--workdir", str(destination.resolve())], check=True)
        return
    destination.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f"apmcm-{protocol}-") as temp:
        work = Path(temp) / "experiment"
        work.mkdir()
        lag = protocol == "delay_identification"
        code = work / "code" if lag else work
        inputs = work / ("inputs" if lag else "input")
        code.mkdir(exist_ok=True)
        inputs.mkdir()
        for file in sorted((ROOT / settings["curve_input"]).glob("curves_*.csv")):
            shutil.copy2(file, inputs / file.name)
        if len(list(inputs.glob("curves_*.csv"))) != 10:
            raise ValueError("Expected ten canonical process curve inputs")
        if lag:
            spectra = ROOT / settings["spectra_input"]
            shutil.copy2(spectra, inputs / spectra.name)
        for script in settings["scripts"]:
            shutil.copy2(source / script, code / script)
            subprocess.run([sys.executable, "-X", "utf8", str(code / script)], cwd=work, check=True)
        if lag:
            for folder in ("tables", "figures"):
                shutil.copytree(work / "rerun_output" / folder, destination / folder, dirs_exist_ok=True)
            summarize_delay(destination / "tables")
            for file in (work / "rerun_output").glob("*.md"):
                shutil.copy2(file, destination / file.name)
        else:
            for extension, folder in (("csv", "tables"), ("png", "figures")):
                target = destination / folder
                target.mkdir(exist_ok=True)
                for file in work.glob(f"*.{extension}"):
                    shutil.copy2(file, target / file.name)

def summarize_delay(tables):
    """Publish the compact table used by the final delay report."""
    full = pd.read_csv(tables / "full_series_lag_diagnostics.csv")
    support = pd.read_csv(tables / "90s_window_support_counts.csv").set_index(["scenario", "gas"])
    risks = pd.read_csv(tables / "90s_window_loho_models.csv").pivot_table(
        index=["scenario", "gas"], columns="model", values="MAE_s", aggfunc="mean"
    )
    rows = []
    for key, g in full.groupby(["scenario", "gas"], sort=False):
        window, risk = support.loc[key], risks.loc[key]
        rows.append({
            "scenario": key[0],
            "gas": key[1],
            "full_lag_median_s": g.best_abs_lag_s.median(),
            "full_peak_abs_r_median": g.best_abs_r.abs().median(),
            "full_BH_q_le_0_1_heats": int((g.BH_q <= .1).sum()),
            "bootstrap_ci_width_median_s": (g.bootstrap_lag_ci_high_s - g.bootstrap_lag_ci_low_s).median(),
            "90s_windows": int(window.windows),
            "90s_threshold_informative": int(window.threshold_informative),
            "90s_BH_supported": int(window.BH_supported),
            "LOHO_fixed_MAE_s": risk.fixed_median,
            "LOHO_time_MAE_s": risk.time_only,
            "LOHO_flow_MAE_s": risk.flow_only,
            "LOHO_time_flow_MAE_s": risk.time_plus_flow,
            "LOHO_time_flow_noninc_MAE_s": risk.time_plus_flow_noninc,
        })
    pd.DataFrame(rows).to_csv(tables / "final_delay_evidence_summary.csv", index=False)
