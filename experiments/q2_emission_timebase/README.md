# Q2.2 cumulative-emission timebase sensitivity

From the repository root, with the corrected official process data extracted and dependencies installed:

    python experiments/q2_emission_timebase/run.py

Optional --input and --output arguments specify the process CSV and destination. By default the source is data/interim/extracted/process.csv and output is results/q2/emission_timebase/rerun/. The input must include row,O,CO2,CO,attachment,heat; 10 official heats and 13,589 records. H0 uses 0.5 s for all rows, H1 uses 1 s only for attachment 5 (6C02057). Original process data are not modified.

The original experiment's fitting and scoring logic is unchanged; only the input/output interface was adapted. On a process CSV reconstructed from the original experiment's ten NPZ arrays, the adapted entry reproduced summary.csv, folds.csv, and predictions.csv byte for byte. The repository's full process CSV was not independently rerun in this commit.

See [report](../../reports/working/Q2_2_EMISSION_TIMEBASE_REPORT.md) and [frozen scores](../../results/q2/emission_timebase/).
