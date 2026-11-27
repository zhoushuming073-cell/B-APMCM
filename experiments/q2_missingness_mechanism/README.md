# Q2: Fifth-heat 6C02057 missingness-mechanism experiment

## Research question

Heat 6C02057 has 775 rows, while the other five heats in the same 6C family have 1,343–1,555 rows. Its cumulative oxygen increase is about 3,536.60 m³. At the official 0.5 s per row, its mean oxygen flow is about 9.14 m³/s; at 1.0 s per row it is about 4.57 m³/s, closer to the other heats. The experiment distinguishes:

- H0: 0.5 s per row, complete records, genuinely higher oxygen flow;
- H1: approximately 1.0 s per row in the fifth heat, with complete records;
- H2: native 0.5 s acquisition with irregular row loss in the final table;
- sensor-specific loss: spectrum, flue-gas, and oxygen streams independently lose rows and are then compressed/re-paired by row order.

H2 is a hypothesis to test, not a preset conclusion.

## Design

1. Use only the six 6C heats to avoid mixing the distinct 6B and 6C wavelength grids.
2. Treat 6C02057 as the target and the other five 6C heats as control sequences.
3. Apply regular systematic thinning, random thinning, and bursty Markov thinning with mean cycle lengths 4, 8, 16, and 32 rows. The primary comparison retains exactly 775 rows in each control heat.
4. Run a retention-rate sensitivity using the rate implied by the median oxygen flow of the five control heats. This estimate depends on the assumption that the fifth heat's true mean oxygen flow is comparable; it is not independent evidence.
5. Compare CO₂, CO, cumulative-oxygen update events, scalar spectral features, and full-spectrum adjacent-frame similarity. Separately compare cross-sensor dependence under a shared global mask and independently compressed sensor streams.
6. Summarize Monte Carlo replicates within each control heat before aggregating across the five control heats. Replicates are not counted as additional independent heats.

## Run

From the repository root:

```bash
python experiments/q2_missingness_mechanism/run.py \\
  --data-dir data/raw/official/workbooks \\
  --out results/q2/missingness_mechanism \\
  --n-rep 100 --cross-rep 100 --shape-rep 20 --seed 20261128
```

Run the additional spectral–gas alignment diagnostic:

```bash
python experiments/q2_missingness_mechanism/spectral_gas_alignment.py \\
  --spectra-dir data/interim/extracted/spectra \\
  --process-csv data/interim/extracted/process.csv \\
  --out results/q2/missingness_mechanism/spectral_gas_alignment \\
  --rep 25 --max-lag 60 --seed 20261128
```

For faster reruns of the main experiment when combined per-heat NPZ caches already exist, add `--cache-dir PATH`; the official XLSB files can still be supplied with `--data-dir` so their hashes are recorded. Requires Python 3.12+, NumPy 2.x, pandas 2.x, pyxlsb 1.x, and SciPy 1.12+ for the spectral–gas diagnostic. Raw inputs are read-only; all outputs go to separate directories.

## Published results

The compact, inspectable result set is in `results/q2/missingness_mechanism/`, including `REPORT.md`, the data audit, sensor-family match summary, full-spectrum shape comparison, and cross-sensor alignment summaries. The full Monte Carlo draws can be regenerated with the commands above and are not required to reproduce the summary tables.

## Interpretation limits

- Five control heats are a small reference sample. Monte Carlo repetitions describe outcomes conditional on those five sequences; they do not create a larger real-world sample.
- Heat-to-heat process heterogeneity can mimic sampling differences.
- Without independent device timestamps, row indices cannot recover actual acquisition times. Cross-sensor dependence can support or weaken a row-alignment hypothesis, but cannot prove physical synchronization.
- This experiment does not automatically replace the official H0 timeline with H2. H0 and H1 remain sensitivity scenarios until independent timestamps or stronger physical constraints become available.
