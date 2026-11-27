# 第一问封存输出

本目录保存十炉最终复算的轻量级冻结摘要。完整 NPZ 情景、逐点饱和审计和图形文件由复算命令生成至 Git 忽略的 `results/q1/rerun/finalization/`，避免在 Git 历史中重复存储大体积矩阵。

- `q1_full_reconstruction_summary.csv`：十炉筛查、异常修复和饱和处理总账。
- `q1_plot_manifest.csv`：十炉三维图的时间/波长显示分辨率及上限标记数量。
- `q1_saturation_structure.csv`：饱和观测的相邻帧可用性和连通区域统计。
- `q1_saturation_simulation_medians.csv`：人工截顶对照的炉次合并中位误差。

复算和核验命令见[实验入口说明](../../../experiments/q1_spectral_reconstruction/README.md)与[第一问封存记录](../../../reports/final/q1/第一问封存记录.md)。
