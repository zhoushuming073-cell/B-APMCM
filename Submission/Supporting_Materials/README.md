# Supporting Materials — Q1 & Q2 (staging snapshot)

**Source commit:** `0fe80e406617e5881bf78efd6ba059ef2fb13302`. Q1 and Q2 are frozen on 2026-11-28. Q3 is intentionally excluded until formally frozen. This is a staging area, not yet the final ZIP or a completed submission.

## Model and evidence index

| Part | Frozen model and entry point | Evidence |
|---|---|---|
| Q1 | Multiscale screening, conservative two-sided repair, conditional hard-clipping lower bound, bilinear continuous reconstruction: `experiments/q1_spectral_reconstruction/run.py` | `reports/final/q1/第一问封存记录.md`, `results/q1/finalization/` |
| Q2.1 | Nonnegative integral cubic B-spline; 2-second local means + PCHIP and analytic integration: `experiments/q2_finalization/run.py` (calls V4 oxygen flow and V1 anomaly screening) | `reports/final/Q2_1_FINAL_FREEZE.md`, `results/q2/final_freeze_2026-11-28/`, validation tables |
| Q2.2 | Oxygen-dose relationships and stage-dependent marginal emissions: `experiments/q2_relationship_study/run.py` | `reports/final/Q2_2_FINAL_RELATIONSHIP_STUDY.md`, `results/q2/relationship_study_final/` |
| Q2.3 | Apparent correlation peak, flow groups, robustness and identifiability checks: relationship-study, lag-robustness and stage-matched-lag scripts | `reports/final/Q2_3_FINAL_VALIDATION_AND_FREEZE.md`, `reports/final/Q2_3_MECHANISM_LITERATURE_AND_CLOSURE.md`, `results/q2/q2_3_final_validation_2026-11-28/` |

H0 is the primary timebase; H1 is a sensitivity scenario only. Q2.3 reports an apparent response lag, not a verified physical transport delay. No monotonic decrease with flow is established. The +8/+10 s CO2 matching discrepancy is unresolved and should not be quoted as a definitive number.

## Reproduction (from this directory, in a clean local copy)

The layout intentionally preserves the source repository's relative paths to avoid changing model code. Install the original Windows/Python 3.13 lock (see requirements.lock.txt), then install this local package. Use the official corrected competition Data.zip, supplied separately by the organizers. Original data and extracted raw-equivalent caches are not part of the submission.

```bash
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python scripts/prepare_data.py /path/to/official/Data.zip
python scripts/verify_data.py
python -m pytest tests/test_q1.py tests/test_data_loader.py -q
python experiments/q1_spectral_reconstruction/run.py
python scripts/verify_q1_outputs.py
python experiments/q2_finalization/run.py
python experiments/q2_relationship_study/run.py
python experiments/q2_lag_robustness/run.py
python experiments/q2_stage_matched_lag/v1/run.py
python experiments/q2_stage_matched_lag/v2/power_identifiability.py
```

Run in a disposable working copy because some original Q2 scripts write into frozen-result directory names. Preserve the committed tables as references. The data-preparation script is new packaging glue, not a modified research model; its end-to-end reproduction and CSV byte equality have not yet been independently verified.

## Outstanding before final ZIP

- Recover the separate `Q2_3_Final_Validation_and_Freeze.zip` artifact, especially `run_q23_final_validation.py` and its supporting scripts, environment and full simulation outputs. The current source tree only contains its README, final report and compact CSV summaries. Do not treat those as a complete reproducibility package.
- Add the actual final Q1 hard-clip 10-heat figures, generated NPZ/processed outputs if needed and within the size budget (do not include raw-equivalent caches or original observations solely as copies), and final Q2 figures/curves consistent with the frozen models. The newer Q1 figures are git-ignored rerun products, not the older tracked figures. Filter regenerated material to comply with original-data exclusion.
- Verify regenerated results against archived summaries and record any differences. No fresh experiment run has been performed by this staging operation.
- Add Q3 after freeze; reconcile final paper figure/table references and AI disclosure separately.
- Final supporting file: `B{TeamNumber}fj.zip`, up to 200 MB, optional on the competition upload page. Do not include raw Data.zip, XLSB workbooks, raw spectral caches, Git LFS pointers, .git, virtual environments or superseded research archive. Check anonymity and remove identifying metadata and source repository links from the final submitted copy where appropriate. Actual submission must remain private during the competition.

`SOURCE_MANIFEST.tsv` records each verbatim copied source blob, path and byte count. Some copied historical robustness reports link to omitted source artifacts; the authoritative model/evidence index above and final README should be used for navigation. No original model or result files have been edited.
