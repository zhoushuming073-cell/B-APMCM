# Q2.3 Stage-matched lag experiment — Round 3

This completed post-freeze supplement tests whether higher oxygen-lance flow is associated with shorter CO/CO2 effective response lag after matching windows within the same heat on cumulative oxygen progress. It does not overwrite or supersede the frozen Q2.2/Q2.3 analyses.

## Reproduce

From the repository root, ensure the frozen Q2.1 curves exist. If needed:

\`\`\`bash
python experiments/q2_finalization/run.py
\`\`\`

Then run:

\`\`\`bash
python experiments/q2_stage_matched_lag/v1/run.py
\`\`\`

Optional explicit paths:

\`\`\`bash
python experiments/q2_stage_matched_lag/v1/run.py \\
  --input results/q2/final_freeze_2026-11-28/final_curves.csv \\
  --output results/q2/stage_matched_lag_2026-11-28
\`\`\`

Dependencies: numpy, pandas, scipy.

## Outputs

The script writes:
- \`lag_windows.csv\`: window-level Pearson and Spearman lag estimates.
- \`matched_pairs.csv\`: one-to-one high/low flow window matches.
- \`heat_level_effects.csv\`: heat-level matched effects.
- \`matching_summary.csv\`: cluster-bootstrap intervals, exact sign tests and sensitivity checks.
- \`manifest.json\`: input fingerprint and experiment settings.

A committed snapshot of the complete outputs and input-reconstruction QA is stored under \`results/q2/stage_matched_lag_2026-11-28/\`.

## Matching and inference

Within each heat and gas, low-flow windows are strictly below the within-heat 25th percentile and high-flow windows strictly above the 75th percentile. Quantile ties are excluded. High- and low-flow windows are paired one-to-one without replacement by minimum absolute cumulative-oxygen-progress distance. The primary diagnostic caliper is 0.15 (15 percentage points of full-heat oxygen progress); 0.05, 0.10, 0.20 and alternative quantile cutoffs are sensitivity checks.

The outcome is \`high-flow lag - low-flow lag\`; a negative value supports shorter lag at high flow. Inference is clustered at the heat level. A missing match indicates lack of common support, not evidence of no effect.

## Interpretation

The primary H0 / 90 s / Pearson / 25th–75th percentile comparison produced seven matched pairs across five of ten heats for each gas. The heat-level median lag difference was +4 s for CO and +8 s for CO2; both heat-cluster bootstrap intervals included zero. Only one of five matched heats had a negative effect for either gas. No matches were available at progress calipers of 0.05 or 0.10.

This is an exploratory post-freeze diagnostic. It does not identify a causal effect and does not establish that higher flow lengthens lag. See the round-three report for the full sensitivity results and evidence boundaries.
