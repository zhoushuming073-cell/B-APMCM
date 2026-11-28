# Q3 Low-Frequency Flue-Gas Feedback Audit

**Date:** 2026-11-28  
**Branch:** research/q3-feedback-audit-20261128  
**Scope:** Reproduce the frozen V5 q2_fixed causal spectral model; audit delayed-observation feedback; evaluate analyzer sampling intervals and score-start sensitivity.

## 1. Data and reproduction

The audit uses the official Data.zip supplied with the project.

SHA-256: cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727

The archive contains 10 heats and 13,589 original time points. A compact cache was rebuilt from the raw workbooks using the Q3 feature definitions in the repository: 32 equal wavelength bands with log-mean intensity, 32 line-contrast features, cumulative oxygen, and CO/CO2 observations.

The reconstructed-cache run of the frozen V5 q2_fixed configuration matches the repository's frozen per-heat R² values to five decimal places for the 120 s score-start / 60 s feedback setting. This establishes a reproducible route from the supplied raw archive to the frozen V5 results; it does not imply byte identity with the separate historical NPZ archive referenced by the old README.

Validation protocol:
- Outer validation: leave one heat out (10 outer folds).
- Inner selection: leave one of the remaining nine heats out.
- Lag selection and model selection use training heats only.
- Feedback interval and gain are selected within the inner folds.
- Feedback uses only analyzer observations available no later than the forecast publication time.

The principal score is the macro-average of the 10 per-heat R² values. Pooled R² is reported separately and must not be substituted for macro R².

## 2. Feedback mechanism and causality

For a detector-clock observation at time s, the measured gas value corresponds to the process state delayed by the transport/analysis lag. The model prediction at the same detector timestamp is built from spectral features at s − theta. Their difference is therefore a residual for the same delayed process state.

At publication time t, the correction uses the most recent scheduled analyzer observation s satisfying s ≤ t:

    y_fb(t) = max(0, y_hat(t) + gamma * [c(s) - y_hat(s)])

The V5 index is anchored at the first scored detector timestamp. With 0.5 s sampling, lag theta seconds and feedback interval Delta seconds, the selected historical index is:

    k_s = floor((k - 2*theta)/(2*Delta)) * (2*Delta)

It is used only when it is nonnegative. Thus the selected detector timestamp is no later than the publication time. The implementation does not use future analyzer readings. The first correction becomes available only after the publication clock reaches the first scored detector timestamp.

This is a valid online delayed-observation correction provided the plant actually exposes those past analyzer readings at the stated timestamps. It is not a pure spectral-only prediction; feedback and no-feedback scores must be reported separately.

## 3. Feedback-interval sweep

All results below use the V5 q2_fixed model family and nested leave-one-heat-out selection. The reduction factor is relative to the raw 0.5 s analyzer sampling interval.

### Common score start: 60 s

| Feedback interval | CO macro R² | CO2 macro R² | CO pooled R² | CO2 pooled R² | Sampling reduction |
|---:|---:|---:|---:|---:|---:|
| 5 s | 0.8828 | **0.7797** | 0.9220 | **0.8232** | 10x |
| 15 s | 0.8774 | 0.7645 | **0.9141** | 0.8124 | 30x |
| 30 s | 0.8808 | 0.7327 | 0.9097 | 0.7900 | 60x |
| 60 s | **0.8871** | 0.7183 | 0.9024 | 0.7788 | 120x |

The best precision-only interval is 5 s. However, 15 s preserves most of the performance while reducing analyzer sampling frequency by 30 times. A 30 s interval reduces it by 60 times, with a larger CO2 penalty. The 60 s interval is competitive for CO but weaker for CO2.

### Score-start sensitivity

| Score start | Feedback interval | CO macro R² | CO2 macro R² | CO pooled R² | CO2 pooled R² |
|---:|---:|---:|---:|---:|---:|
| 60 s | 15 s | 0.8774 | 0.7645 | 0.9141 | 0.8124 |
| 90 s | 15 s | **0.8908** | 0.6607 | 0.8967 | 0.7548 |
| 120 s | 15 s | 0.7785 | 0.6278 | 0.8844 | 0.8082 |

The score window has a material effect. It is not legitimate to select a score start solely because it maximizes R²; the window must be justified by analyzer validity and the task's intended prediction horizon.

## 4. Empirical response-onset diagnostic

A descriptive onset diagnostic was calculated as the first time the gas signal remained above 1 m³/s for 10 consecutive seconds:
- CO2: latest onset across the ten heats was approximately 52 s.
- CO: latest onset across the ten heats was approximately 115 s.

Using a 0.5 m³/s threshold gave latest onsets of approximately 49 s for CO2 and 109.5 s for CO. This supports 60 s as a conservative common CO2 response window and 120 s as a conservative CO response window. The threshold is an empirical diagnostic, not a threshold specified by the contest problem.

A gas-specific evaluation window (CO2 from 60 s; CO from 120 s) is defensible if explicitly justified, while a common 60 s window is useful as a full-system sensitivity result. Both should be disclosed rather than silently choosing the more favorable one. Under that gas-specific onset window with 15 s feedback, the macro R² is 0.7645 for CO2 and 0.7785 for CO; the corresponding pooled R² values are 0.8124 and 0.8844. These two gas scores use different time supports and must not be presented as a common-window comparison.

## 5. Other exploratory checks

### Residual update variants

Last-residual hold, rolling median of the last 3/5 residuals, exponential decay, and a shrinkage affine residual calibration were compared under nested validation. No complex update rule produced a consistent macro-R² improvement over the simple last-residual hold. The additional flexibility did not justify the complexity.

### Converter-group expert averaging

Experts were grouped by the heat-ID prefixes 6B and 6C, then compared with the global expert average and mixtures. The grouped expert was not consistently selected by inner validation and did not improve both gases. The prefix-to-converter mapping itself should be confirmed before interpreting this as a definitive test of physical converter grouping.

### Eight-feature proxy inspired by the external method note

A proxy using log brightness, four broad-band color ratios, K/Na line-to-continuum contrasts, trailing flicker, and cumulative-oxygen splines was implemented. It did not outperform V5. This is not an exact reproduction of the external method: the supplied screenshots do not specify its exact wavelength intervals, line-strength definitions, or standardization scope.

### Per-heat prefix normalization

Using the first 60 s of each heat for causal feature centering/scaling did not produce a material overall improvement. The current V5 normalization/feature representation remains preferable on the evidence tested.

## 6. Interpretation and recommendation

1. Historical flue-gas feedback is compatible with a real-time forecasting interpretation when the observation timestamp is no later than the forecast publication time. It must be described as a feedback-enhanced predictor, not a spectral-only predictor.
2. The V5 feedback mechanism is already implemented; the main new finding is the interval and score-window audit.
3. For a common 60 s score start, 15 s feedback is the recommended engineering compromise: CO macro R² 0.8774 and CO2 macro R² 0.7645, with 30x lower analyzer sampling frequency. If accuracy is prioritized over sampling reduction, 5 s is stronger. If CO is prioritized, 60 s remains competitive but harms CO2.
4. The 60/90/120 s score-start sensitivity must be reported. A conservative onset-based evaluation can use CO2 from 60 s and CO from 120 s; do not claim that the 60 s common-window score is directly comparable with a 120 s evaluation.
5. The worst heat remains a significant limitation. At a common 60 s start and 15 s feedback, minimum per-heat R² is 0.4331 for CO and 0.1622 for CO2. The model has not achieved uniformly high per-heat generalization comparable to the external claim.
6. The external report's R² protocol is not yet known to be identical. Its per-heat versus pooled definition, score window, excluded heats, and standardization scope must be obtained before a direct numerical comparison.

## 7. Limitations

- Only ten heats are available; heat-level bootstrap intervals are descriptive and are not corrected for the multiple interval comparisons.
- Pooled R² benefits from between-heat level variation and can be much higher than macro R².
- The response-onset threshold is a transparent diagnostic, not an official contest rule.
- The external method was only approximated from screenshots; its exact feature engineering and validation procedure remain unknown.
- The V5 consolidated repository module previously had not been rerun end-to-end; this audit reconstructed the cache and reproduced the frozen q2_fixed per-heat results, but the historical original NPZ archive was not available for byte-level comparison.

## 8. Files

- experiments/q3_feedback_audit/prepare_cache.py: rebuild compact features from official Data.zip.
- experiments/q3_feedback_audit/run_feedback_sweep.py: nested LOHO feedback and score-start sweep.
- results/q3/feedback_audit/summary.csv: macro, pooled, paired-difference and sampling-reduction summary.
- results/q3/feedback_audit/fold_metrics_60s_15s.csv: per-heat baseline and 15 s feedback scores for the common 60 s window.
