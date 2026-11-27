# Frozen results: Q2.2 cumulative-emission timebase sensitivity

- summary.csv: 16 records covering H0/H1, two gases, two models, in-sample and LOHO; both macro and pooled R².
- folds.csv: 160 individual heat scores. Both files match the original rerun's Git blob hashes exactly.
- predictions.csv: the reproduction entry regenerates 24,080 direct mean-curve comparison points; not duplicated in Git. The original file SHA-256 is 127ce62a6263893d5153664a0a5fc80d2b8cebdc09b061bc27d1ffbbcbb2f7ad.

Direct mean curve: 301 equally spaced common net-oxygen points from 0 to 3306.63 m³; cumulative emission integrated from first record, repeated oxygen levels resolved by keeping the last record. Sigmoid integral: original absolute oxygen coordinate, fixed row index >=120, time integration reset at evaluation start. These two domains differ; their R² values are not strictly comparable. H1 is a hypothesis, not a verified correction to official metadata.

[Reproduction](../../../experiments/q2_emission_timebase/) · [Research report](../../../reports/working/Q2_2_EMISSION_TIMEBASE_REPORT.md)
