# Q1 spectral transition study: candidate spectral markers around CO2-to-CO transition

## 1. Question and scope
This exploratory study asks which spectral bands or candidate atomic lines change consistently near the gas-defined rise in CO share, whether current spectra identify the contemporaneous stage in an unseen heat, and whether current spectral level or short-term rate is more informative. It does not claim future prediction.

## 2. Data and safeguards
- Ten heats, official 0.5 s per row, 2047 spectral channels per row.
- The two instrument wavelength grids were mapped to a common 0.338 nm grid over approximately 340.535–1023.303 nm.
- Raw spectra are the primary analysis. Saturated channels are masked in line-effect estimation; where a classifier requires complete input, clipping is filled by same-row interpolation over wavelength only.
- Q1 reconstruction is retrospective sensitivity analysis only: anomaly repair can use future frames (up to 12 rows) and calibration uses the first 55% of a heat; saturation repair uses adjacent frames.
- Cross-heat inference treats heats, not frames, as the independent units. Classification uses leave-one-heat-out (LOHO), with feature selection fitted within each training fold.

## 3. Gas-defined transition
CO share is CO/(CO+CO2). CO and CO2 are median-filtered, and the share is smoothed before defining the first time it remains at or above 0.5 for 10 s while total gas exceeds 0.5. Across the ten heats, the median transition time is 234.75 s; range 112.5–297.5 s. Spectral contrasts compare [-40,-10] s before and [+10,+40] s after the transition.

## 4. Strongest candidate: K I doublet neighbourhood
The strongest cross-heat candidate is the 765–770 nm band, aligned with the NIST K I resonance lines at 766.490 and 769.896 nm. For the pre-specified wide doublet contrast, raw spectra increased in 9/10 heats at the 50% CO-share transition; mean paired change was 0.01534 log-contrast units, with unadjusted paired t-test p=0.001712. The 40% threshold showed a similar direction (9/10); the 60% threshold was weaker (7/10, p=0.494). This is a candidate spectral marker, not proof of potassium concentration change or a causal chemical mechanism.

The peak-centroid estimates were within about 0.2 nm of the NIST wavelengths, comparable to the common-grid channel spacing (~0.338 nm). Instrument response, temperature, excitation, self-absorption and wavelength calibration remain alternative explanations.

## 5. Current-stage classification (not forecasting)
Two distinct LOHO classification pipelines were run and must not be conflated:
- Pre-specified candidate-line pipeline: raw K_doublet_wide_only mean balanced accuracy 0.619 and mean AUC 0.741.
- Separate full-transition-window comparison pipeline: K_doublet mean balanced accuracy 0.611 and AUC 0.701; 2047-channel full_line mean balanced accuracy 0.647 and AUC 0.716.

These differ in feature construction and sampling pipeline; both assess contemporaneous state near an event-defined transition. They do not establish that K leads the gas response.

For backward 5 s rates, the candidate-line classifier was near chance (raw K_doublet_wide_only balanced accuracy 0.512, AUC 0.519); the separate full-window K_doublet_rate comparison was 0.470 / 0.473. Current level carries more stage information than short-term rate in these tests.

## 6. Other candidate regions
- 775–780 nm: 9/10 heats showed a decrease in the 5 nm local line-contrast band (BH q=0.0249); it overlaps the O I 777 nm triplet, but the narrow O I line-group test was weaker and does not confirm the assignment.
- 600–605 nm: all ten heats decreased; BH q=0.00666, but no reliable atomic-line attribution has been established.
- 475–480 nm: 9/10 heats increased, but BH q=0.0634, above the 5% FDR threshold.
- Na I D region around 589 nm did not pass the multiple-testing threshold.

## 7. Conclusions and boundaries
Supported: a repeatable local spectral-contrast change occurs near the gas-defined CO-share transition; the strongest candidate is the K I doublet neighbourhood; raw and Q1-retrospective results point in a similar direction; current spectral level has modest contemporaneous stage information.

Not supported: potassium abundance change; a causal K-to-CO chemical mechanism; reliable early warning from 5 s spectral rate; future CO-share prediction; or a definitive O I / Na I assignment.

## 8. Next experiment
Use raw spectra only in a strict rolling-origin design. At each forecast origin, use only current and past spectra, past gas observations, instantaneous oxygen-lance flow estimated by backward differencing cumulative oxygen, and cumulative oxygen to predict future CO-share changes at 10, 20 and 30 s. Compare a gas/process baseline against baseline plus K I features, and test whether any gain survives leave-one-heat-out validation.
