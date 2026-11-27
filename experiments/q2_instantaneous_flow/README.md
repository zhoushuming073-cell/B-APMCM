# Q2.2 instantaneous flow V1–V3

Run from repository root after installing the repository dependencies:

    python experiments/q2_instantaneous_flow/v1/run.py --input data/interim/extracted/process.csv --out results/q2/instantaneous_flow/rerun/v1
    python experiments/q2_instantaneous_flow/v2/run.py --input data/interim/extracted/process.csv --out results/q2/instantaneous_flow/rerun/v2
    python experiments/q2_instantaneous_flow/v3/run.py --input data/interim/extracted/process.csv --out results/q2/instantaneous_flow/rerun/v3

All three scripts are byte-identical to their original ZIP versions (Git blob hashes verified). Always specify --input and --out; the preserved original defaults refer to the local experiment layout. The original shared process CSV was constructed from the official ten-heat cache (Git blob 94bc24258ad638bc5ac9c8d6bc544524757a4ec4). The repository process CSV is the intended equivalent input; exact archive reproduction can instead use process.csv in any original ZIP. This commit does not independently rerun the entire repository data pipeline.

Use --out outside the frozen results. See [frozen data](../../results/q2/instantaneous_flow/) and [research report](../../reports/working/Q2_2_INSTANTANEOUS_FLOW_REPORT.md).
