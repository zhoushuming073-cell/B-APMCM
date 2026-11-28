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
- [Complete reproducibility bundle](Q2_3_Stage_Matched_Lag_Round3.zip)
- manifest.json: input fingerprints and settings.
- input_reconstruction_qa.json: source and reconstructed-curve checks.

The ZIP contains the full window-level lag table, matched-pair table, heat-level effects, complete sensitivity summary, script, report, manifest and QA record. The complete output tables can also be regenerated from the frozen Q2.1 curves using the experiment entry point.
