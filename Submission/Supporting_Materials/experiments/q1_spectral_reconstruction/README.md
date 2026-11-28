# 第一问：全域光谱重构与封存复算

本入口重建十炉原生波长网格上的异常筛查、保守修复情景、饱和情景及三维曲面。正式数据为更正版十炉，共 13,589 帧、每帧 2,047 个原生波长通道。

## 主流程

在仓库根目录执行：

```bash
python scripts/verify_data.py --rebuild-cache
python experiments/q1_spectral_reconstruction/run.py
python scripts/verify_q1_outputs.py
python -m pytest tests/test_q1.py -q
```

默认输出至 `results/q1/rerun/finalization/`，不会覆盖 `results/q1/diagnostics/` 中的冻结实验表。也可使用 `--output PATH` 指定其他目录。

指定炉次复算示例：

```bash
python experiments/q1_spectral_reconstruction/run.py --heat 6B03205 --heat 6C02057
```

## 输出解释

每炉 `scenarios/{HEAT}_scenario.npz` 保存：

- `raw`：原始观测，不改写；
- `scenario_anomaly_only`：仅批准异常候选替换后的情景；
- `scenario_template`：加入邻帧饱和模板估值的情景；
- `scenario_hard_clip_lower_bound`：假定 65,535 为硬截顶时施加下界约束的条件情景；
- `screen`、`approved`、`saturated`、`saturation_repaired`：逐点来源掩码；
- `saturation_template_values`、`saturation_template_below_ceiling`、`saturation_template_above_ceiling`、`saturation_unresolved`：饱和模板的逐点估值及可用性。

历史字段 `scenario` 保留为 `scenario_template` 的兼容别名。它不是硬截顶下界情景。不能将邻帧模板估值表述为实测峰顶。

`tables/q1_full_reconstruction_summary.csv` 是十炉总账；`tables/q1_saturation_audit.csv` 逐点记录上限观测、模板估值和条件下界。 `scripts/verify_q1_outputs.py` 会重新检查掩码关系、数据哈希、输出形状及三种情景之间的值来源，并从 NPZ 文件重建这两张汇总表。

## 三维图

主 PDF `figures/all_heats_spectral_surface.pdf` 使用 `scenario_hard_clip_lower_bound`，红色标记记录值达到 65,535 的位置。波长轴默认不降采样（2,047 个通道全部显示）；时间轴每 10 帧取样一次（5 s）以控制 3D 网面渲染成本。节点分析仍在完整时间—波长网格上进行。图是显示抽样，不是重新拟合。

需要展示未施加硬截顶下界的邻帧模板敏感性图时：

```bash
python experiments/q1_spectral_reconstruction/run.py --template-sensitivity
```

单炉绘图只处理命令行指定的炉次，不会误将目录中其他炉次的旧情景混入所选 PDF。

## 独立诊断

以下命令在同一入口下运行峰形/注入诊断、七通道峰遮挡指标复核和饱和结构/人工截顶模拟：

```bash
python experiments/q1_spectral_reconstruction/run.py --diagnostics
```

其中 `peak_study.py` 和 `peak_fit_metrics.py` 使用既有诊断协议；`saturation_study.py` 已从历史 `/mnt/data` 绝对路径脚本迁移为项目内模块，输入统一来自 `data/interim/extracted/spectra/`，输出写入 `results/q1/rerun/saturation_study/`。

## 证据边界

- 人工注入召回率不是自然异常识别准确率；原始数据没有人工异常标签。
- 65,535 的仪器截顶机制尚无设备文档独立确认。
- 邻帧模板是条件估值；真实截顶前峰顶没有独立真值。
- 双向邻帧使用属于离线回顾性重构，不是实时预测流程。
- 连续曲面按原生节点做双线性张量积插值，节点之间的数值不是新增观测。

详见[第一问封存记录](../../reports/final/q1/第一问封存记录.md)。
