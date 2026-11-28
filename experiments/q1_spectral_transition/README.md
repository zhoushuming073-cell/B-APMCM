# Q1 spectral transition study (10 heats)

## Purpose
Identify spectral bands and candidate atomic lines associated with the gas-defined CO2-to-CO share transition. Compare current spectral level with strictly backward spectral rates.

## Reproduction
1. Extract the original Data.zip so that the ten Attachment N. Heat No. ...xlsb workbooks are in Data/.
2. Install dependencies from requirements.txt.
3. Run source/run_all.py.
4. Review results/REPORT.md, results/tables/, and results/figures/.

The Q1 reconstruction runs one heat per process to limit memory. Intermediate arrays are stored in cache/. Raw spectra remain the primary source for any real-time or causal claim. Q1 reconstructed spectra are retrospective sensitivity data because repair and calibration use future frames or full-heat information.

## Validation safeguards
- Heat is the independent unit for cross-heat inference.
- Classification is leave-one-heat-out; feature selection is performed within training folds.
- Raw spectral features use the current frame only; temporal rates use backward differences.
- Event-window classification describes current-stage discrimination near a gas-defined transition, not future prediction.
- Two wavelength-axis groups are aligned to a common wavelength grid before pooled comparisons.

## Main outputs
See results/REPORT.md and results/tables/. The next experiment should use raw spectra only and a strict rolling-origin design.
