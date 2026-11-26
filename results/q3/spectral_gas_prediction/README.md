# 第三问冻结实验结果（2026-11-27）

19 模型族 × 两气体 × 两时基 × 10 炉 = 760 项留炉评分。854,924 条逐点预测全部复算验证。这里保留可直接查看的冻结文本表，不重复上传原始数据或大体积结果。

- macro_metrics.csv：76 个模型族/气体/情景汇总，十炉等权 MAE、RMSE、R²、最差炉。
- key_fold_mae.csv：峰形光谱集成与同结构纯过程集成的十炉 MAE（六位小数）。
- paired_effects.csv：主要消融的配对炉均改善、改善炉数、描述性重采样区间。
- lags.csv：40 个外层折的训练侧提前量。
- lag_negative_control.csv：全炉事后描述性循环移位对照（不参与训练）。

完整逐时点预测、全部逐炉指标/超参数、图表由 [实验入口](../../../experiments/q3_spectral_gas_prediction/) 生成到 rerun/。详见 [模型与实验报告](../../../reports/working/Q3_SPECTRAL_GAS_MODEL_REPORT.md)。冻结结果为探索性实验，不代表现场验收。
