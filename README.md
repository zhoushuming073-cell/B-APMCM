# 2026 APMCM B 题研究仓库

转炉火焰光谱连续重构、供氧与烟气响应时滞分析。正式分析使用更正版十个炉次，共 13,589 帧。本仓库保存数据审计、模型代码、实验协议、冻结结果和研究报告。

## 当前研究边界

最新研究状态、推荐模型、证据等级、待决问题及全部专项索引见 **[RESEARCH_STATUS.md](RESEARCH_STATUS.md)**（2026-11-28 更新）。

- **第一问：已封存。** 十炉全量复算，保守异常修复、饱和情景和双线性连续重构已冻结；模板值不是真实峰顶。
- **第二问 2.1：已封存。** 非负积分 B 样条估计有效供氧流量，k2/PCHIP 重构烟气并解析积分；51 个异常候选仅标记，原始数据不修改。
- **第二问关系研究：已封存。** 使用 Q2.1 冻结曲线完成累计耗氧—烟气流量/累计排放关系及流量分组时滞研究；不做预测。高流量对应更短时滞的方向未获稳健支持。
- **第二问时滞：已并入最终关系研究。** 90/180 秒窗口的有效时滞估计存在尺度敏感性；结果不支持强制拟合单调递减时滞函数。
- **第三问：探索性。** 已完成光谱融合、曲线迁移、残差树、神经网络、梯度提升和物理过程实验；CO 梯度提升与 CO₂ 残差树为当前候选，非最终模型。光谱增量须按相同信息集及评分域审查。
- **第五炉：未决。** H0 仍为主情景，H1 为敏感性；H2 不规则缺行获得描述性支持，真实秒级时间轴及跨传感器缺帧同步性未识别。

## 快速开始

以下命令均在仓库根目录运行。依赖锁按 Windows / Python 3.13 环境解析。

```powershell
git lfs pull
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python scripts/verify_data.py
python -m pytest
```

已封存 Q1 和 Q2.1 可分别复算：

```powershell
python experiments/q1_spectral_reconstruction/run.py
python scripts/verify_q1_outputs.py
python experiments/q2_finalization/run.py
python experiments/q2_relationship_study/run.py
```

查看既有批量实验入口（仅覆盖脚本注册的专项，不含全部新增实验）：

```powershell
python scripts/run_all.py --help
python scripts/run_all.py q2_optimized
python scripts/run_all.py q3_spectral_gas_prediction
python scripts/run_all.py
python scripts/build_figures.py
```

批量复算输出到各实验的 `rerun/` 目录，重绘图形输出到 `results/rebuilt_figures/`。**尚未独立复核的探索性实验不会自动加入默认批量运行**；Q2.3 预白化实验如需单独检查，可运行：

```powershell
python experiments/q2_prewhitened_delay/run.py
```

该脚本默认输出到 `results/q2/prewhitened_delay/rerun/`，不会覆盖已归档的 `tables/` 结果。具体边界见 [实验说明](experiments/q2_prewhitened_delay/README.md) 和 [Q2.3 研究报告](reports/working/Q2_3_PREWHITENED_DELAY_REPORT.md)。

## 研究导航

| 研究线 | 报告 / 入口 | 冻结结果 |
| --- | --- | --- |
| 数据与审计 | [正式数据](data/raw/official/Data.zip) · [审计表](data/metadata/audit.csv) · [SHA-256](data/metadata/checksums.sha256) | [数据审计说明](archive/original_full_experiment_package/APMCM_Full_Experiments/README.md) |
| 第一问：光谱重构 | [封存记录](reports/final/q1/第一问封存记录.md) · [复算入口](experiments/q1_spectral_reconstruction/README.md) | [十炉复算总账](results/q1/finalization/q1_full_reconstruction_summary.csv) · [窄峰遮挡指标](results/q1/tables/q1_peak_mask_fit_metrics.csv) |
| 第二问 2.1：正式封存 | [封存报告](reports/final/Q2_1_FINAL_FREEZE.md) · [复算入口](experiments/q2_finalization/) | [冻结 QA](results/q2/final_freeze_2026-11-28/) · [重构验证](results/q2/validation_2026-11-28/) |
| 第二问关系与时滞研究：最终封存 | [最终研究报告](reports/final/Q2_2_FINAL_RELATIONSHIP_STUDY.md) · [复算入口](experiments/q2_relationship_study/) | [紧凑结果](results/q2/relationship_study_final/) |
| 第二问 2.1：解析供氧 V4、异常筛查 | [流量 V4](experiments/q2_instantaneous_flow/v4/) · [异常检测 V1](experiments/q2_anomaly_review/v1/) · [更新报告](reports/working/Q2_2_FLOW_ANOMALY_UPDATE.md) | [V4 结果](results/q2/instantaneous_flow/v4/) · [异常候选](results/q2/anomaly_review/v1/) |
| 第二问：第五炉缺失机制 | [实验说明](experiments/q2_missingness_mechanism/) · [研究报告](results/q2/missingness_mechanism/REPORT.md) | [模拟结果](results/q2/missingness_mechanism/) |
| 第二问：对齐比较 | [实验入口](experiments/q2_aligned_comparison/) | [结果](results/q2/aligned_comparison/) |
| 第二问：优化实验 | [优化实验报告](reports/final/q2_optimized/Q2_report.pdf) | [90 秒窗口证据](results/q2/delay_identification/tables/final_delay_evidence_summary.csv) |
| 第二问：时基审计 | [时基审计报告](reports/final/timebase_audit/APMCM_Timebase_Audit_Report.pdf) | [时基审计结果](results/q2/timebase_audit/) |
| 第二问 2.2：累计排放时基（历史敏感性） | [历史报告](reports/working/Q2_2_EMISSION_TIMEBASE_REPORT.md) · [实验入口](experiments/q2_emission_timebase/) | [历史结果](results/q2/emission_timebase/) |
| 第二问 2.2：瞬时流量 V1–V3（历史模型探索） | [历史报告](reports/working/Q2_2_INSTANTANEOUS_FLOW_REPORT.md) · [实验入口](experiments/q2_instantaneous_flow/) | [历史结果](results/q2/instantaneous_flow/) |
| 第二问 2.3：原时滞识别 | [实验入口](experiments/q2_delay_identification/) | [结果](results/q2/delay_identification/) |
| 第二问 2.3：预白化时滞（探索性） | [研究报告](reports/working/Q2_3_PREWHITENED_DELAY_REPORT.md) · [实验说明](experiments/q2_prewhitened_delay/README.md) | [结果表](results/q2/prewhitened_delay/tables/) |
| 第三问：光谱—烟气预测 | [模型报告](reports/working/Q3_SPECTRAL_GAS_MODEL_REPORT.md) · [实验入口](experiments/q3_spectral_gas_prediction/) | [冻结结果](results/q3/spectral_gas_prediction/) |
| 第三问：B 炉专项 V1–V3（探索性） | [模型报告](reports/working/Q3_B_HOLDOUT_MODEL_REPORT.md) · [实验入口](experiments/q3_b_holdout/) | [冻结评分](results/q3/b_holdout/) |
| 第三问：V5 光谱与低频烟气反馈（探索性） | [研究报告](reports/working/Q3_V5_CAUSAL_SPECTRAL_REPORT.md) · [实验入口](experiments/q3_v5_causal_spectral/) | [冻结评分](results/q3/v5_causal_spectral/) |
| 第三问：曲线迁移、残差树与神经网络 | [方差报告](reports/working/第三问异质性曲线机器学习与方差.md) · [运行说明](experiments/q3_curve_transfer/README.md) | [逐炉指标](results/q3/curve_transfer/folds.csv) · [方差汇总](results/q3/curve_transfer/variance.csv) |
| 第三问：梯度提升追加实验 | [20 个 R² 报告](reports/working/第三问梯度提升追加实验20个R2.md) · [运行说明](experiments/q3_boosted_curve/README.md) | [20 个 R²](results/q3/boosted_curve/twenty_r2.csv) · [20 张图](results/q3/boosted_curve/figures/README.md) |
| 第三问：十炉物理过程 | [研究报告](reports/working/q3_physical_process/REPORT.md) · [实验说明](experiments/q3_physical_process/) | [结果](results/q3/physical_process/) |

最终交付：[论文](reports/submission/paper/) · [摘要页](reports/submission/summary_sheet/) · [控制页](reports/submission/control_sheet/)。目前为待整合交付区域，已有正式封存材料在 [reports/final](reports/final/)。

研究过程文档：[实验进度](reports/working/实验进度.md) · [数学模型与实验证据](reports/working/2026%20APMCM%20B题——第三版数学模型与实验证据.md) · [第三问项目指南](reports/working/THIRD_QUESTION_PROJECT_GUIDE.md) · [第三问研究记录](reports/working/第三问研究.docx)。

## 目录结构

| 目录 | 职责 |
| --- | --- |
| `configs/` | 数据位置、模型参数、实验输入与复算配置 |
| `docs/` | 赛题资料、官方规范、方法与验证协议 |
| `data/` | 原始附件、提取数据、处理数据及审计 |
| `src/apmcm/` | 数据读取和 Q1–Q3 核心代码 |
| `experiments/` | 各专项的独立运行入口与说明 |
| `results/` | 实验输出、冻结表格和图形 |
| `reports/` | 工作报告、最终报告和提交材料 |
| `archive/` | 迁移前交付包、被替代实验及历史材料 |

## Git 协作约定

- `main` 是唯一长期维护和集成分支。
- 专项实验确有隔离需要时，才建立短期分支；不要为每次参数调整长期保留分支。
- 实验完成后通过 Pull Request 合并到 `main`，确认合并后删除临时分支。
- 冻结结果、探索性结果和已验证结论必须明确区分；不覆盖历史结果来制造连续性。
- 新增专项时同步更新 [研究状态总账](RESEARCH_STATUS.md) 及本页导航，注明实验输入、代码、结果、证据边界和复跑状态。
