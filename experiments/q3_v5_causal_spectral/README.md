# V5 causal spectral LOHO

Install the repository dependencies. Extract the original Q3_V5_Causal_Spectral_LOHO.zip provided in the research conversation. Use its exact q3_eda NPZ cache.

From repository root:

    python src/apmcm/q3/v5_causal_spectral.py --cache /path/to/q3_eda --out results/q3/v5_causal_spectral/rerun_spectral --lag-scheme spectral_lag
    python src/apmcm/q3/v5_causal_spectral.py --cache /path/to/q3_eda --out results/q3/v5_causal_spectral/rerun_q2_fixed --lag-scheme q2_fixed

The repository module consolidates the two original scripts by a lag-scheme CLI option; the model logic and settings are retained. The consolidated module has not been rerun end to end during this commit. Frozen results are unchanged. Original complete archive SHA-256: 8d96954291cdec7c507336053e8ec0a19cbd378a67b44b88bfb8d298ae2d3010. Full-precision per-fold scores, inner-fold tables, predictions and cache remain in that archive, not duplicated in Git. See [report](../../reports/working/Q3_V5_CAUSAL_SPECTRAL_REPORT.md).
