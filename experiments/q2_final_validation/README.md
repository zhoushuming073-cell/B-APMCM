# Q2.3 Final Validation

Final Q2.3 validation is complete. The detailed report and compact CSV outputs are stored under:
- reports/final/Q2_3_FINAL_VALIDATION_AND_FREEZE.md
- results/q2/q2_3_final_validation_2026-11-28/

The full reproducibility package contains the one-command entry point, frozen Q2.1 flow profile, scripts, full 100-replicate semi-synthetic outputs, figures, and environment pins. Command:

```bash
python run_q23_final_validation.py --data-zip Data.zip --out q23_final_output --simulations 100 --seed 20261128
```

The output directory must be empty. The frozen Q2.1 profile is reused rather than rerunning its expensive cross-validation search; all 20 heat/scenario flow summaries match the official main-branch reference. The 100-replicate simulation is split into deterministic batches of ten.

Final status: **带明确限制封存**. No fifth Q2.3 model-search round is planned.
