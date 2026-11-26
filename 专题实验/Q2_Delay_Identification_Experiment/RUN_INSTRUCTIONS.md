# Reproducing the delay experiments

The package includes the ten per-heat curve CSVs in `inputs/`. The main spectral diagnostic table is already included. The original `harmonized_spectra.npz` is not bundled because it is about 80 MB; to rerun the spectral-brightness component, place that file in `inputs/`.

Required Python packages: `numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-learn`, and `statsmodels`.

Recommended rerun order from the unzipped package root:

1. `python code/run_q2_delay_experiment.py` — full-series lag profiles, block bootstrap, circular-shift diagnostics, and baseline local windows.
2. `python code/postprocess_q2_delay.py` — summarize full-series results and create baseline figures.
3. `python code/run_q2_time_adjusted.py` then `python code/postprocess_q2_time_adjusted.py` — exploratory 120 s non-overlap sensitivity. This is retained for audit only; it is not the final inference because several heats have too few usable windows.
4. `python code/run_q2_window_robustness.py` — final 90 s non-overlap window tests, window-level circular-shift controls, BH adjustment, and leave-one-heat-out comparisons.

Rerun outputs are written to `rerun_output/` so the delivered result set is not overwritten. The scripts use bundle-relative input paths.
