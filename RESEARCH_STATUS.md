# 研究状态总账 | RESEARCH_STATUS

更新日期：2026-11-28（北京时间）。本次研究证据快照：[4d23a08](https://github.com/zhoushuming073-cell/B-APMCM/commit/4d23a08b52d260e4ea158f5f8ac2cc3206a37fcf)。本文是研究决策导航，不替代各专项报告、原始数据或冻结结果；后续有新实验时同步更新本页及 [README](README.md)。

## 1. 当前决策与证据等级

- **F（已封存）**：已有明确冻结决定、模型及归档验证结果；不意味着物理真值已独立确认。
- **E（探索性实证）**：有实际实验或归档评分，但使用反复探索的十炉，尚非独立确认；复跑状态以原报告为准。
- **H（待验证假设）**：有机制线索或先验，但不可视为已证实。
- **P（待完成）**：尚需决策、复核或形成交付物。

| 研究线 | 等级 | 当前决策 | 首要证据 |
| --- | --- | --- | --- |
| Q1 光谱连续重构 | F | 主算法、参数、情景及绘图口径已封存 | [封存记录](reports/final/q1/第一问封存记录.md) |
| Q2.1 六类过程曲线 | F | 非负积分 B 样条 + k2/PCHIP；原始烟气仅标记异常，不自动修改 | [封存报告](reports/final/Q2_1_FINAL_FREEZE.md) |
| Q2.2 耗氧与排放/瞬时流量关系 | E | 简洁共享曲线及 Sigmoid 积分为基准；复杂异质性模型未稳定取胜，待收敛 | [累计量报告](reports/working/Q2_2_EMISSION_TIMEBASE_REPORT.md) · [瞬时流量报告](reports/working/Q2_2_INSTANTANEOUS_FLOW_REPORT.md) |
| Q2.3 有效响应时滞 | E/H | 保留固定时滞基准；流量条件变时滞尚不稳健 | [预白化报告](reports/working/Q2_3_PREWHITENED_DELAY_REPORT.md) |
| Q3 实时烟气预测 | E | 优先比较 CO 梯度提升、CO₂ 融合残差树，保留共享曲线及纯过程对照；未封存 | [梯度提升报告](reports/working/第三问梯度提升追加实验20个R2.md) · [曲线迁移报告](reports/working/第三问异质性曲线机器学习与方差.md) |
| 第五炉时基与缺失 | E/H | H0 主情景、H1 敏感性；H2 不规则缺失为优先研究假设，真实时间轴未识别 | [缺失机制报告](results/q2/missingness_mechanism/REPORT.md) |
| 竞赛论文与提交包 | P | 优先将已封存结果和未决边界整合成论文 | [论文目录](reports/submission/paper/) · [摘要页](reports/submission/summary_sheet/) · [控制页](reports/submission/control_sheet/) |

## 2. 数据与共同约束

正式数据：[Data.zip](data/raw/official/Data.zip)，SHA-256 为 `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`。10 炉、13,589 帧/行，每帧 2,047 个原生波长通道。数据以[审计表](data/metadata/audit.csv)和[校验值](data/metadata/checksums.sha256)为准。

- H0：题面 0.5 s/行，全部十炉；H1：仅第五炉 6C02057 按 1 s/行，其他炉仍为 0.5 s/行。H2：不规则缺行，真实时间戳不可识别。不得改写原始数据或把 H1/H2 当成正式时间轴。
- 第五炉 775 行，净耗氧 3536.60 m³；在对照工况可比假设下，估计保留率约 52.3%，非实测丢帧率。随机和短簇状缺失比规则降采样更匹配部分烟气和氧计数器统计；跨光谱、烟气的缺失同步性未识别。
- 独立样本单位是炉次，不是时点；同一批十炉的多轮模型选择不构成独立验证。Q3 发布时刻只能使用此前已可获得的信息；离线双向重构、最终耗氧、整炉归一化及未来数据不能直接用于实时预测。
- R² 区分炉均与合并，累计量与瞬时量、H0/H1、评分域、预测提前量和反馈条件不得混比。

## 3. Q1：已封存

- 多尺度时空残差筛查、双向模板保守修复、原生波长双线性连续曲面；保留 raw、仅异常修复、模板饱和估值、硬截顶下界情景。
- 十炉全量复算：51,236 个候选、5,074 个批准修复点；14,805 个上限读数，其中 14,373 获模板估值、432 未解决。合成筛查点召回 88.24%，非真实标注准确率。
- 正式图采用硬截顶下界情景，65,535 上限点作标记；模板估值不等于真实峰顶，双向修复不等于在线算法。后续仅做论文表述和排版，非必要不重新搜参。
- [冻结结果及曲面](results/q1/finalization/) · [入口](experiments/q1_spectral_reconstruction/) · [诊断](results/q1/diagnostics/) · [封存记录](reports/final/q1/第一问封存记录.md)。

## 4. Q2：过程曲线、关系与时滞

### Q2.1 已封存（六条曲线）
非负三次 B 样条直接表示有效瞬时供氧，解析积分拟合累计氧；CO₂/CO 采用 2 秒局部均值 + PCHIP（k2），解析积分得到累计排放。所有原始烟气观测保留，51 个候选（3 单点、10 多点、38 阶跃）仅标记待审。

| 指标 | H0 | H1 |
| --- | ---: | ---: |
| 累计氧阻断验证 MAE (m³) | 2.141 | 2.182 |
| CO₂ k2 MAE (m³/s) | 0.04708 | 0.04711 |
| CO k2 MAE (m³/s) | 0.12829 | 0.13189 |

20/20 拟合成功、非负、累计单调、质量平衡检查通过（最大数值差 0.0773 m³）；累计排放与原始梯形积分差异不超过 0.021%。无独立高频流量真值、无自然异常标签；不得把合成注入召回当作真实准确率。原始数据零修改。

[正式冻结报告](reports/final/Q2_1_FINAL_FREEZE.md) · [复算入口](experiments/q2_finalization/) · [冻结 QA 与汇总](results/q2/final_freeze_2026-11-28/) · [重构验证](results/q2/validation_2026-11-28/) · [解析流量 V4](experiments/q2_instantaneous_flow/v4/) · [异常检测](experiments/q2_anomaly_review/v1/)。

### Q2.2 尚待最终决策
- 累计量：直接平均曲线 H0 LOHO R²（CO₂/CO）0.4926/0.7631，H1 为 0.9662/0.9402；Sigmoid 后积分为 0.9645/0.9539（两时基逐炉分数不变，评分域不同，不直接判优）。高累计量 R² 不证明瞬时关系准确或 H1 成立。
- 瞬时量：H1 共享曲线 CO₂ 0.6596、CO 0.8457，V1 CO 0.8557；静态/动态异质性未稳定优于简单基准。H1 数值不能冒充 H0 主结果。
- 建议：将冻结的 Q2.1 曲线与旧实验的原始标签协议明确区分；统一 H0/H1、评价区间、误差及排放终点指标后选择论文模型。联合碳流/氧化比例可作探索，不硬设未知物理比例。
- [累计量研究](reports/working/Q2_2_EMISSION_TIMEBASE_REPORT.md) · [瞬时量 V1–V3](reports/working/Q2_2_INSTANTANEOUS_FLOW_REPORT.md) · [流量/异常更新](reports/working/Q2_2_FLOW_ANOMALY_UPDATE.md) · [累计结果](results/q2/emission_timebase/) · [瞬时结果](results/q2/instantaneous_flow/)。

### Q2.3 尚待最终决策
90/180 秒窗口、去趋势、AR(1) 预白化及 Pearson/Spearman/距离相关平方结果存在尺度与时基敏感性。H0 90 秒窗口无 BH q≤0.1 的诊断窗口；180 秒仅部分 Pearson CO₂ 有支持，其他指标未复现。H0 下 12 项 LOHO 线性流量模型均劣于固定均值基准。相关峰值包含反应、输运、仪器响应，不是独立物理运输时间真值。

保留固定时滞基准、变时滞物理先验与局限，避免挑选单一有利窗口。预白化 Python 移植未独立端到端复跑；原结果表为归档证据。

[既有优化报告](reports/final/q2_optimized/Q2_report.pdf) · [对齐比较](results/q2/aligned_comparison/) · [时滞证据](results/q2/delay_identification/) · [预白化报告](reports/working/Q2_3_PREWHITENED_DELAY_REPORT.md) · [窗口结果](results/q2/prewhitened_delay/tables/) · [文献笔记](docs/methodology/q2_delay_optimization_literature.md)。

## 5. Q3：模型候选与淘汰依据（尚未封存）

下表按 CO₂ / CO 列示。Q3 主候选实验均为同一批十炉上的探索，跨轮数值须核查评分域和提前量；以下不构成独立确认。

| 实验 | 结果/定位 | 证据 |
| --- | --- | --- |
| 初始 19 模型族双时基嵌套留炉 | 760 项评分；光谱增量未稳定证实，保留作历史基准 | [报告](reports/working/Q3_SPECTRAL_GAS_MODEL_REPORT.md) |
| B 炉 V1–V3 | 单炉专项；V2 B 炉 R² 0.7574/0.9032，不代表十炉泛化 | [报告](reports/working/Q3_B_HOLDOUT_MODEL_REPORT.md) |
| V5 因果光谱 | Q2 固定时滞改造：无反馈融合 0.4633/0.6776；60 s 反馈 0.5790/0.7760；反馈是不同信息集，单列情景 | [报告](reports/working/Q3_V5_CAUSAL_SPECTRAL_REPORT.md) |
| 曲线迁移 | H0 共享曲线 0.4719/0.7155，融合残差树 0.4862/0.7209，融合网络 0.4930/0.7168；树跨三情景更稳健 | [报告](reports/working/第三问异质性曲线机器学习与方差.md) |
| 梯度提升 | H0 融合 0.4519/0.7308；纯过程同结构 0.3841/0.6991；CO 候选优先，CO₂ 保留残差树；第五炉 CO₂ 为负 | [报告](reports/working/第三问梯度提升追加实验20个R2.md) |
| 十炉物理探索 | CO₂ 起升先于 CO：10/10 炉，阈值定义下中位领先 77.5 s；光谱单独识别阶段未稳定超过时间基准；非因果结论 | [报告](reports/working/q3_physical_process/REPORT.md) |

建议采用共享耗氧进度曲线 + 因果光谱残差校正为论文候选架构。严格分别报告纯过程、仅光谱（若有相同协议结果）、融合和可用反馈情景；优先针对 CO₂ 跨炉不稳定及第五/第九炉，检查模型选择是否依赖外层成绩。无独立新炉时明确限制，不宣称已证明光谱普遍增益。物理探索阶段标签、居中平滑及整炉统计不得泄漏进实时模型。

[曲线迁移结果](results/q3/curve_transfer/) · [梯度提升 20 个 R²](results/q3/boosted_curve/twenty_r2.csv) · [20 张图](results/q3/boosted_curve/figures/README.md) · [物理过程结果](results/q3/physical_process/)。

## 6. 实验与结果完整索引

以下覆盖当前快照全部一级专项实验目录；版本子目录按功能索引。入口存在不代表当前环境已复跑，冻结大表可能由脚本生成或保留于历史包。

| 专项 | 实验入口/说明 | 结果 |
| --- | --- | --- |
| Q1 光谱重构 | [入口](experiments/q1_spectral_reconstruction/) | [Q1](results/q1/) |
| Q2 对齐比较 | [入口](experiments/q2_aligned_comparison/) | [结果](results/q2/aligned_comparison/) |
| Q2 异常筛查 | [V1](experiments/q2_anomaly_review/v1/) | [结果](results/q2/anomaly_review/) |
| Q2 时滞识别 | [入口](experiments/q2_delay_identification/) | [结果](results/q2/delay_identification/) |
| Q2 排放时基 | [入口](experiments/q2_emission_timebase/) | [结果](results/q2/emission_timebase/) |
| Q2.1 正式封存 | [入口](experiments/q2_finalization/) | [结果](results/q2/final_freeze_2026-11-28/) |
| Q2 瞬时流量/异质性 | [V1–V3](experiments/q2_instantaneous_flow/) · [V4（Q2.1）](experiments/q2_instantaneous_flow/v4/) | [结果](results/q2/instantaneous_flow/) |
| Q2 缺失机制 | [入口](experiments/q2_missingness_mechanism/) | [结果](results/q2/missingness_mechanism/) |
| Q2 优化实验 | [入口](experiments/q2_optimized/) | [结果](results/q2/optimized/) |
| Q2 预白化时滞 | [入口](experiments/q2_prewhitened_delay/) | [结果](results/q2/prewhitened_delay/) |
| Q2 时基审计 | [入口](experiments/q2_timebase_audit/) | [结果](results/q2/timebase_audit/) |
| Q3 B 炉 | [入口](experiments/q3_b_holdout/) | [结果](results/q3/b_holdout/) |
| Q3 梯度提升 | [入口](experiments/q3_boosted_curve/) | [结果](results/q3/boosted_curve/) |
| Q3 曲线迁移 | [入口](experiments/q3_curve_transfer/) | [结果](results/q3/curve_transfer/) |
| Q3 物理过程 | [说明](experiments/q3_physical_process/) | [结果](results/q3/physical_process/) |
| Q3 初始光谱气体预测 | [入口](experiments/q3_spectral_gas_prediction/) | [结果](results/q3/spectral_gas_prediction/) |
| Q3 V5 | [入口](experiments/q3_v5_causal_spectral/) | [结果](results/q3/v5_causal_spectral/) |

## 7. 下一步（建议优先级）

1. P0：锁定 Q2.2、Q2.3 的论文结论和基准；对照 Q2.1 冻结口径，明确 H0 主情景、H1/H2 局限。
2. P0：统一 Q3 公平比较条件，决定每种气体的最终候选与退回模型，报告逐炉分数、误差和失败案例。
3. P0：撰写竞赛论文、摘要与控制页，引用冻结表和明确的探索性结果；不将报告目录中的历史版本当作最新交付。
4. P1：若可获得设备日志/时间戳，核验第五炉；否则保留不可识别结论及 H0/H1/H2 敏感性。不人为重建缺失原始观测。
5. P1：核验关键复跑脚本、图表和最终提交材料；尤其 Q2.3 Python 移植、Q3 已迁移入口与历史结果一致性。

## 8. 后续提交约定

- 实验新增或变更时，记录问题编号、数据版本、信息集、验证协议、运行入口、输出路径、关键结果、证据限制和当前决策；新增专项同步更新第 6 节及首页导航。
- 已封存模型如需改变，新增版本并记录原因，不覆盖历史表格、不把探索性分数写成正式结论。全文报告优先于本总账的摘要数值。
- 结果大文件按既有 LFS/本地复算约定管理；历史资料留在 archive。默认批量脚本仅包含其显式注册的实验，不代表全部专项。
- 本次是文档和索引审查，未独立重新执行全部实验。
