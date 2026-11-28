# Q3 feedback fine search (2026-11-28)

Exploratory, not frozen for submission. This supersedes the invalid 60/90 s early-window claims in the earlier feedback audit (full 60-second prefix unavailable at some publication times).

The two Python source files are byte-identical to the local experiment package (Git blob IDs 9454c758988b7195a21ac2877c00cc078b83861d and 1a5ac7a7be78658f72f313ae06825baf07bf7b0d).

From repository root, build the cache using the previous prepare_cache.py:

```bash
python experiments/q3_feedback_audit/prepare_cache.py --archive data/raw/official/Data.zip --out data/interim/q3_feedback_cache
python experiments/q3_feedback_fine_search/run_fine_search.py --cache data/interim/q3_feedback_cache --out results/q3/feedback_fine_search/rerun120 --score-start-s 120 --prefix-s 60 --workers 4
python experiments/q3_feedback_fine_search/run_fine_search.py --cache data/interim/q3_feedback_cache --out results/q3/feedback_fine_search/rerun60 --score-start-s 60 --prefix-s 0 --workers 4
python experiments/q3_feedback_fine_search/run_expanded_search.py --cache data/interim/q3_feedback_cache --out results/q3/feedback_fine_search/rerun_expanded60 --score-start-s 60 --prefix-s 0 --workers 4
```

Requires NumPy, pandas, SciPy, scikit-learn, threadpoolctl, pyxlsb and repository dependencies. Run on Linux (process pool inherits loaded cache). The archived scripts write macro scores and per-fold results; pooled R2 in the local package was calculated separately after runs. The original archive has 7,560 outer score rows and cache files; this commit preserves exact scripts plus selected rounded numerical tables, not the entire archive. Original archive SHA256: d2011c516c3760313deabef2fdeb1ee097c32ab52a2116bf02c958e4ca0436de. Data.zip SHA256: cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727.

## Method and audit

- Ten outer leave-one-heat-out folds, nine inner folds per outer fold. Equal-heat scores, train-only scaling, train-only lead selection. First frame only in corrected 60 s detector-start experiments; 120 s detector start permits the full 60 s prefix.
- Process template indexed by estimated oxygen progress; initial spectrum predicts template shape/scale; ongoing spectral and process features fit ridge/PLS/heat-expert candidates. Residual feedback uses latest available detector observation, not future target.
- Detector-clock target at t+d, publication at t, feedback observation s <= t. The lead d is selected in training folds from 6,12,...,60 seconds. Centered smoothing is applied only to historical training targets when constructing templates; runtime spectral features use trailing windows.
- Exhaustive 30 intervals (0.5–15 s), coarse gain 0–1.5 step .25, fine gain step .05 among top three coarse base candidates per interval. Expanded alpha grid .01,.03,.1,.3,1,3,10; residual strengths .25,.5,.75,1,1.25. Per-fold inner selection, not globally fixed parameters.
- Interval curves are outer sensitivity comparisons; selecting a favorite after seeing outer scores is exploratory selection. Nested interval rows are the appropriate evaluation of the interval-selection procedure. No independent confirmation cohort.
- Lead estimated with a Q2-inspired training-only criterion, not a direct reuse of the final Q2 variable-lag results. Q1 offline continuous reconstruction is not fully integrated. H0 only; fifth heat unresolved; no hardware validation of actual analyzer duty reduction. The 60 s and 120 s evaluation supports are different. Historical gas feedback is permitted as a plausible interpretation of dual-modal fusion, not an explicit official ruling.
- No evidence from the reviewed indexing of target-future access in corrected experiments; an independent streaming/prefix-invariance test and full end-to-end rerun are still needed for final acceptance.

## Key results

See results/q3/feedback_fine_search. At common detector start 60 s and expanded grid, fusion macro R2 CO/CO2: no feedback .835938/.641512, 0.5 s .894707/.794038, 5 s .893289/.776372. Pooled at 5 s .929786/.822997. The CO fixed-interval curve actually peaks at 14 s (.895621) in expanded60; this is an outer sensitivity maximum, not a validated optimum. CO2 prefix-only prior + feedback at 0.5 s .812420 beats ongoing spectral fusion .794038. At 5 s the gas-only persistence control is .478633/.188279; process+feedback .880239/.761884. Thus feedback and process templates explain most accuracy, with modest incremental ongoing spectral contribution for the compared families. The fifth heat has .454266/.200625 at 5 s. Do not claim universal high R2 or independent confirmation.
