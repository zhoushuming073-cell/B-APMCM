# Q2.3 Stage-matched lag — Round 3 results

Completed post-freeze supplementary diagnostic. This result set does not replace the frozen Q2.2/Q2.3 outputs.

## Primary result

H0, 90 s, Pearson association, within-heat 25th/75th flow tails, oxygen-progress caliper 0.15:

| Gas | Matched pairs | Heats | Heat-level median lag difference | 95% heat-cluster bootstrap interval | Heats with shorter lag at high flow |
| --- | ---: | ---: | ---: | ---: | ---: |
| CO | 7 | 5/10 | +4 s | [−4, +6] s | 1/5 |
| CO2 | 7 | 5/10 | +8 s | [−26, +20] s | 1/5 |

The effect is high-flow lag minus low-flow lag. Negative values would support the hypothesized shorter lag at higher flow. The intervals include zero.

## Main limitation

No matches were available at oxygen-progress calipers 0.05 or 0.10. At the primary 0.15 caliper, only five of ten heats contributed matches, and the median progress gap was 0.134. This is relaxed stage matching, not an exact same-stage comparison. The result is not causal evidence and does not establish that high flow lengthens lag.

## Files

- [Round-three report](../../../reports/final/Q2_3_STAGE_MATCHED_LAG_ROUND3_ADDENDUM.md)
- [Experiment code and run instructions](../../../experiments/q2_stage_matched_lag/v1/)
- key_results.csv: compact H0/H1 primary and sensitivity results.
- manifest.json: input fingerprints and settings.
- input_reconstruction_qa.json: source and reconstructed-curve checks.

The full per-window lag estimates, matched pairs and heat-level effects are generated locally by the experiment script and are not committed as large derived tables. This follows the repository convention for reproducible outputs: retain compact evidence and regenerate detailed tables from the frozen Q2.1 curves when needed.

The rebuilt curve hash differs from the repository-recorded frozen CSV hash, while the numerical headline summaries were reproduced. This is numerical-level verification, not a byte-identical reproduction; see input_reconstruction_qa.json and the report.
