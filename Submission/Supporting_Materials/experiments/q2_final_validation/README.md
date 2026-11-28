# Q2.3 — Frozen final-validation reproduction (submission edition)

The scientific model and final-validation protocol are unchanged. All eight original Python files, dependency pins and 25 retained archive members are preserved. The excluded member, `inputs/process.csv`, is a direct raw-equivalent extraction of the competition workbooks; the original runner generates that input from official Data.zip at runtime.

From the Supporting_Materials root (Python 3.13):

```bash
python experiments/q2_final_validation/restore_package.py --out q23_reproduction
python -m pip install -r q23_reproduction/requirements.txt
python q23_reproduction/run_q23_final_validation.py --data-zip /path/to/official/Data.zip --out q23_final_output --simulations 100 --seed 20261128
```

Both output directories must be empty. The submission archive consists of ten parts and is 3322278 bytes; SHA-256: `ff402e5482ec70d58cfd677c14af8a9e9ed669c59d19216b01cd80c5dec08117`. The upstream archive SHA-256 is `a72c7e7a596ee20f7fb016c9384e1c51ad92b40613a76200dedc5d342c5d2c0d`. The new manifest retains the original member hashes for all 25 members; the unchanged restore code validates the whole archive and each member. The retained frozen reference curves contain model-derived values alongside observation columns needed for historical numerical audits. Do not distribute locally generated raw input or cache files in the final submission ZIP.

See [original package instructions](PACKAGE_README.md), [final report](../../reports/final/Q2_3_FINAL_VALIDATION_AND_FREEZE.md), [mechanism and closure](../../reports/final/Q2_3_MECHANISM_LITERATURE_AND_CLOSURE.md), and [matching audit](../../results/q2/q2_3_final_validation_2026-11-28/MATCHING_SUMMARY_CONSISTENCY_AUDIT.md). Use heat-level median CO +4 s and CO2 +8 s; pair-level CO2 +10 s is descriptive only. No refit or new model search was conducted. This session verified the upstream archive SHA-256 and retained ZIP entries, not a fresh full restoration or 100-replicate execution.
