# APMCM 2026 Problem B — Full 10-Heat Experimental Report

**Data version:** user-uploaded `Data.zip`  
**ZIP SHA-256:** `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`  
**Scope:** Questions 1 and 2. Question 3 prediction model was not part of this run.  
**Evidence label:** exploratory experiments; synthetic injection and masked-observation tests are not ground-truth anomaly labels.

## 1. Data audit

The archive contains 10 XLSB workbooks, one worksheet per heat, with 2,047 spectral channels and three process columns (cumulative oxygen, CO₂, CO). There are 13,589 total time points and 27,816,683 spectral observations. All sequence indices are continuous from 1 to the final row within each heat. No missing spectral or process cells were found. Sampling interval is 0.5 s per the problem statement.

| Heat | Rows | Duration (s) | Saturated 65,535 points | Oxygen flat steps | CO₂ zero rows | CO zero rows |
|---|---:|---:|---:|---:|---:|---:|
| 6B03205 | 1,479 | 739.0 | 0 | 779 | 0 | 95 |
| 6B03206 | 1,311 | 655.0 | 13 | 633 | 0 | 92 |
| 6B03270 | 1,413 | 706.0 | 0 | 725 | 77 | 95 |
| 6B03271 | 1,420 | 709.5 | 10 | 714 | 43 | 101 |
| 6C02057 | 775 | 387.0 | 4,168 | 242 | 0 | 32 |
| 6C02058 | 1,343 | 671.0 | 5,805 | 661 | 0 | 95 |
| 6C02059 | 1,555 | 777.0 | 100 | 796 | 0 | 111 |
| 6C02060 | 1,373 | 686.0 | 326 | 689 | 26 | 84 |
| 6C02107 | 1,459 | 729.0 | 2,110 | 755 | 0 | 30 |
| 6C02108 | 1,461 | 730.0 | 2,273 | 741 | 0 | 98 |

There are two exact wavelength grids: the four 6B heats share one grid (340.535–1027.963 nm); the six 6C heats share another (339.390–1023.303 nm). The within-family grids are exactly identical. Their common support is 340.535–1023.109 nm (2,030 points on the 6B reference grid). Linear round-trip interpolation from the 6C grid to the common grid and back has MAE 28.6–40.2 raw intensity units on sampled rows. Thus, native grids should be retained for per-family modeling and plotting; cross-family pooling requires explicit harmonization and must not equate column indices with wavelength.

There are 14,805 exact 65,535 readings. Every saturated pixel has at least one horizontally adjacent saturated wavelength channel. Across saturated pixels, 68.08% have both adjacent time frames unsaturated and 31.92% have one adjacent frame unsaturated; none lack both adjacent frames. This supports temporal reference repair as a conditional scenario, not proof of the instrument's clipping mechanism.

The cumulative oxygen series has 6,735 flat adjacent steps and zero negative steps. Median positive increments are 4.72–4.75 m³ across heats. The median update gap is 1.0 s for nine heats and 0.5 s for 6C02057. This is consistent with a sample-and-hold / refresh process, but does not prove the exact meter mechanism. CO₂ has 146 zero rows and CO has 833; zeros were retained as observations.

## 2. Question 1 experiments

### 2.1 Multi-scale anomaly screening

For each heat, the first 55% of frames calibrate per-channel robust scales and a 99.9th-percentile score threshold; a 40-frame guard separates calibration and late test segments. Each test segment contains 144 synthetic anomaly events (1,104 altered points per heat; 1,440 events and 11,040 points in total across ten heats). Events include positive/negative perturbations, one- and multi-channel spectral widths, and one-/multi-frame durations at 5%, 15%, and 40% amplitude.

| Detector | Mean point recall | Mean event recall |
|---|---:|---:|
| Temporal Hampel | 2.12% | 5.42% |
| Prior soft joint residual | 24.30% | 58.06% |
| Multi-scale local residual | **88.24%** | **79.86%** |

For the multi-scale detector, mean event recall by injected amplitude is 51.7% (5%), 88.3% (15%), and 99.6% (40%). It is the best candidate-screening method among those tested. However, its actual late-segment candidate rate averages about 0.42%, above the 0.1% calibration budget, indicating temporal distribution drift. The detector must not be used as an automatic deletion rule.

### 2.2 Conservative two-sided repair gate

The screen is followed by two local shoulder-corrected templates (previous and next frame). A point is repaired only when both templates are available, agree on the correction direction, the minimum departure exceeds 0.12 log-intensity units, and the template disagreement is less than half the departure (equivalently departure > 2 × template gap). The global gate was selected on injected calibration segments, then evaluated on late test segments.

Across ten heats, test-set total MAE fell from 3.065 to 2.217 raw intensity units on average; the test injected-point MAE fell from 1,713 to 694; mean point repair recall was 47.8% and event repair recall 61.6%. Total MAE improved in 7/10 heats. This is evidence that conservative two-sided gating can improve synthetic corruption recovery, but there is still nonzero collateral modification on unaltered points and no independent human anomaly labels. Therefore, this gate produces a **provisional repair scenario**, not a confirmed cleaned truth.

Applying the global gate to all raw data produced 51,236 screened candidates and 5,074 approved anomaly-repair points. Normal observations were retained. Full per-heat masks and repair scenarios are stored in `cleaned_scenarios/`.

### 2.3 Cross-heat peak reconstruction and calibration groups

A leave-one-heat-out masked-peak experiment evaluated 598 algorithmically selected local-maximum candidates per masking width. The held-out row's peak region was hidden; basis coefficients were fitted only from visible wavelengths. Rank-12 basis trained on other heats in the same calibration family outperformed a pooled basis trained on all other nine heats in 8/10, 9/10, and 9/10 heats for widths 3, 7, and 11 channels respectively.

| Mask width | Pooled all-9 rank-12 MAE | Same-family rank-12 MAE | Same-family rank-18 native-grid MAE |
|---|---:|---:|---:|
| 3 channels | 191.1 | 163.1 | 157.9 |
| 7 channels | 175.9 | 150.4 | 145.9 |
| 11 channels | 165.9 | 145.0 | 139.4 |

The rank-18 result is exploratory and was not selected in a separate nested hyperparameter loop; rank 12 remains the more conservative reference. The main finding is that **family-specific bases are preferable to forcing the two wavelength calibration groups into one basis**. Peak candidates are algorithmically selected, not manually verified physical emission lines.

### 2.4 Saturation and censored-signal simulation

Among actual 65,535 points, 14,373/14,805 were assigned a two-sided or one-sided temporal-template estimate; 432 remain unresolved in the provisional scenario. A separate synthetic censoring experiment used 20 unsaturated peak candidates per heat, clipped to 50% and 75% of peak height. Neighbor templates had median peak-height relative error 19.6% at both caps, compared with 50.0% / 25.0% for leaving clipped values unchanged and 50.2% / 31.1% for the tested censored-Gaussian fit. For peak area, the neighbor template was not uniformly best (especially at the 75% cap). This supports conditional peak-height recovery when adjacent frames are unsaturated, but not a universally accurate peak-area reconstruction. The actual clipping threshold and instrument response still require confirmation.

### 2.5 Continuous surfaces

A provisional scenario was generated for all ten heats. It retains the original raw matrices and separately stores candidate masks, approved anomaly repairs, saturation masks, and saturation-template values. Ten 3D surfaces are downsampled for display only; analysis remains on the native grid. First-order tensor-product interpolation is the current continuity-preserving visualization baseline; it does not claim to infer unobserved fast events between samples.

## 3. Question 2 experiments

### 3.1 Oxygen-flow reconstruction

A nonnegative integral B-spline represents latent flow as a nonnegative B-spline expansion and fits the cumulative oxygen record through its integral. Knot spacing and second-difference penalty were selected in inner blocked validation; outer blocks were kept separate for scoring.

| Model | Mean outer-block cumulative-oxygen MAE | Mean outer-block RMSE |
|---|---:|---:|
| Linear interpolation | 2.361 m³ | 2.859 m³ |
| PCHIP | 5.278 m³ | 6.137 m³ |
| Nonnegative integral spline | **2.141 m³** | **2.603 m³** |

The integral spline beats linear interpolation in 7/10 heats, but the gain is modest. There is no independent instantaneous-flow ground truth, so this validates cumulative-signal interpolation, not the physical accuracy of the derivative. In a synthetic sample-and-hold simulation using each heat's empirical update schedule, update-interval average slopes were best when the assumed update schedule was correct; this is conditional on a mechanism that has not been independently confirmed.

The centered 10-second slope gives a useful descriptive flow-support summary. Nine heats have a narrow overlapping P10–P90 support region of approximately 3.89–5.21 m³/s. Heat 6C02057 is distinct, with P10/P50/P90 approximately 5.36/9.07/13.03 m³/s. High-flow support is therefore largely supplied by one heat, creating a serious heat-identity confound for flow-dependent delay estimation.

### 3.2 Gas-signal reconstruction and cumulative emissions

Using the same outer time blocks for all methods, linear interpolation was the strongest simple baseline:

| Gas | Linear MAE | PCHIP MAE | Nonnegative smoothing-spline MAE |
|---|---:|---:|---:|
| CO₂ | **0.0490** | 0.0504 | 0.1034 |
| CO | **0.1441** | 0.1869 | 0.2917 |

PCHIP is a useful shape-preserving sensitivity check; the tested global smoothing spline oversmooths important local variation. Cumulative emissions were computed by trapezoidal integration over the recorded interval; they do not include emissions outside the recording window.

### 3.3 Cumulative oxygen, gas composition, and time-third descriptions

Across heats, the median Spearman correlation of cumulative oxygen with instantaneous CO is 0.817 and with total CO₂+CO is 0.783, while the median correlation with CO₂ alone is −0.155. After first differencing, median correlations are much smaller (CO 0.210, total 0.216, CO₂ −0.036). This demonstrates that raw level correlations are substantially affected by common process trends and should not be interpreted causally.

For descriptive time thirds (not independently identified metallurgical stages), median CO₂ flow falls from 17.06 to 15.75 to 13.07 m³/s, while median CO rises from 3.82 to 40.88 to 43.12 m³/s. Median CO₂ share falls from 0.791 to 0.268 to 0.239. These are descriptive patterns, not a validated stage mechanism.

### 3.4 Delay identification

For each heat and gas, cross-correlation was computed between smoothed oxygen-flow changes and gas-flow changes over 0–40 s. A circular-shift negative control was used; Benjamini–Hochberg correction was applied across 20 heat×gas tests. Six tests had adjusted q < 0.05, but best lags ranged from 0 to 40 s, median absolute correlation was only about 0.123, and several maxima were at the search boundary. More importantly, window-level Spearman correlation between flow and estimated lag was approximately −0.02 for both gases, not supporting the expected monotonic decrease in delay with higher flow. Flow-bin medians were also non-monotone.

A nested leave-one-heat-out prediction comparison used the same 5-second issue time, target, and information set. For CO₂, gas persistence alone achieved MAE 0.311; adding current or historical oxygen-flow features worsened the macro MAE (0.344–0.352). For CO, adding process progress to persistence improved MAE from 0.761 to 0.734 in 9/10 held-out heats; current flow gave 0.736, while historical-flow and lag-grid models were worse (0.780 and 0.767). This is a small predictive increment for CO, not evidence of a physical transport-delay function.

**Decision:** do not publish a single fitted τ(q) curve as established. Report per-heat exploratory lag profiles, boundary hits, flow support, and “insufficient identification” where appropriate. The problem's expected inverse flow-delay relation is not empirically established by this dataset under the tested methods.

## 4. Recommended current model choices

1. Keep original native wavelength grids for all within-family processing and plots.
2. Use multi-scale time–wavelength residuals only to screen candidates; retain a three-tier normal / uncertain / approved-repair mask.
3. For local peak reconstruction, use a family-specific shared basis as a candidate baseline, with rank 12 as the conservative reference and rank 18 as a sensitivity result; compare against two-sided local templates.
4. Treat 65,535 as a possible censoring threshold, retain lower-bound masks, and only use adjacent-frame repair as a conditional scenario. Do not claim the inferred peak top is measured truth.
5. For oxygen flow, use nonnegative integral spline as the primary smooth estimate with linear cumulative interpolation as a benchmark; show uncertainty because there is no flow ground truth.
6. For gas signals, use linear interpolation as the main baseline and PCHIP as a sensitivity check.
7. Do not claim a universal flow-dependent delay. Current evidence supports only exploratory per-heat associations, not a stable τ(q).

## 5. Unresolved items

- No independent human-labelled true anomalies or emission-line catalogue.
- Exact 65,535 instrument saturation mechanism and detector specifications unconfirmed.
- Family and calibration grid are perfectly confounded with heat ID; optical calibration and heat effects cannot be separately identified.
- Instantaneous oxygen-flow ground truth and device timestamps are absent.
- Time delay is not stable under the tested cross-correlation and LOHO prediction designs; higher-flow support is concentrated in one heat.
- Third-question real-time spectral-to-gas prediction was not run in this experiment batch.
- Rank selection above 12 and fully nested family-specific basis selection remain to be completed.

## 6. Reproducibility and output files

See `README.md` for the script inventory, output tables, figure directory, data hashes, and execution order. Original input workbooks were not modified. All repaired spectra are separate scenario files in `cleaned_scenarios/`.
