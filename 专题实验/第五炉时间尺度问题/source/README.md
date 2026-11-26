# APMCM 2026 B — 6C02057 time-base sensitivity audit

## Reproduction

Python 3.10+ is recommended.

```bash
pip install -r requirements.txt
python run_all.py --input Data.zip --workdir apmcm_timebase_results
```

`Data.zip` is read-only. The extractor writes separate XLSB copies under `raw_data/` and numeric NPZ caches under `cache/`. Analysis outputs are written under `results/`; original workbooks are never modified.

If you already have ten XLSB workbooks extracted:

```bash
python run_all.py --input ./Data --workdir apmcm_timebase_results
```

The workflow is split into four scripts so each stage can be rerun independently:

1. `extract_raw_cache.py`: reads each XLSB in an isolated worker process and records file hashes.
2. `analyze_audit.py`: raw sequence/counter audit, spectrum and gas dynamics, oxygen-flow reconstructions, integral-spline CV, cross-correlation delays, and time-scale invariance tables.
3. `fit_lag_models.py`: fixed / flow-conditioned / non-increasing lag models and matched-row leave-one-heat-out finite distributed-lag prediction.
4. `make_figures.py`: generates the 18 comparison figures.

## Hypotheses and interpretation

- H0: 6C02057 uses the official 0.5-s interval.
- H1: 6C02057 uses 1.0 s per row; sensitivity scenario only.
- All other heats remain at 0.5 s in both scenarios.

A counter update gap is counted in rows between positive cumulative-oxygen changes. It is **not** a sampling interval. A row-indexed sequence identifies cumulative volume per row, approximately `q × Δt`; it does not separately identify physical flow `q` and sample interval `Δt` without an external clock or validated physical anchor.

The near-quantum statistic uses the pooled median positive increment as a descriptive reference. It does not establish the instrument's actual pulse resolution. The lag cross-correlations are exploratory and are not ground-truth transport delays. The finite distributed-lag model is predictive, not causal identification.

## Main outputs

- `tables/oxygen_counter_audit.csv`
- `tables/counter_quantization_diagnostics.csv`
- `tables/oxygen_flow_summary.csv`
- `tables/fifth_oxygen_spline_cv.csv`
- `tables/gas_dynamics_summary.csv`
- `tables/gas_integrals_and_balance.csv`
- `tables/spectral_dynamics_summary.csv`
- `tables/lag_estimates.csv`
- `tables/flow_conditioned_lag_models.csv`
- `tables/distributed_lag_LOHO_fair_macro.csv`
- `figures/*.png`

## Validation protocol

The distributed-lag comparison keeps the same raw-row test mask under H0 and H1 (first 83 rows are excluded for lag warm-up), and uses equal total training weight per heat. The official 0.5-s setting remains the primary specification unless independent acquisition timestamps or a separately validated physical constraint supports H1.
