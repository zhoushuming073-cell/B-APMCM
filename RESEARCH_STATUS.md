# 研究状态总账 | RESEARCH_STATUS

更新日期：2026-11-28（北京时间）。本页与 main 分支当前研究状态同步；以 [main 最新提交](https://github.com/zhoushuming073-cell/B-APMCM/commits/main) 为准。本文是研究决策导航，不替代各专项报告、原始数据或冻结结果；后续有新实验时同步更新本页及 [README](README.md)。

## 1. 当前决策与证据等级

- **F（已封存）**：已有明确冻结决定、模型及归档验证结果；不意味着物理真值已独立确认。
- **E（探索性实证）**：有实际实验或归档评分，但使用反复探索的十炉，尚非独立确认；复跑状态以原报告为准。
- **H（待验证假设）**：有机制线索或先验，但不可视为已证实。
- **P（待完成）**：尚需决策、复核或形成交付物。

| 研究线 | 等级 | 当前决策 | 首要证据 |
| --- | --- | --- | --- |
| Q1 光谱连续重构 | F | 主算法、参数、情景及绘图口径已封存 | [封存记录](reports/final/q1/第一问封存记录.md) |
| Q2.1 六类过程曲线 | F | 非负积分 B 样条 + k2/PCHIP；原始烟气仅标记异常，不自动修改 | [封存报告](reports/final/Q2_1_FINAL_FREEZE.md) |
| Q2.2 耗氧—烟气关系与时滞研究 | F/H | 2026-11-28 正式封存；第二轮稳健性检验未验证高流量短时滞；封存后第三轮匹配诊断仍未发现稳健负向证据，且共同支持不足，不重开冻结研究 | [最终研究报告](reports/final/Q2_2_FINAL_RELATIONSHIP_STUDY.md) · [第二轮稳健性检验](reports/final/Q2_2_LAG_HYPOTHESIS_ROBUSTNESS.md) · [第三轮 addendum](reports/final/Q2_3_STAGE_MATCHED_LAG_ROUND3_ADDENDUM.md) |
| Q2.3 有效响应时滞（并入第二问最终研究） | F/H + E（封存后补充） | 原时滞研究随 Q2.2 封存；第三轮同炉同耗氧进度匹配仅 5/10 炉形成主分析匹配，CO/CO₂ 炉次中位时滞差为 +4/+8 s，区间跨零；共同支持不足，不能作因果解释 | [第二轮稳健性检验](reports/final/Q2_2_LAG_HYPOTHESIS_ROBUSTNESS.md) · [第三轮 addendum](reports/final/Q2_3_STAGE_MATCHED_LAG_ROUND3_ADDENDUM.md) |
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

### Q2.2 已完成最终关系研究并正式封存（不做预测）

本轮使用 Q2.1 冻结的十炉六类过程曲线，研究累计耗氧—CO/CO₂ 瞬时流量、累计排放、阶段边际排放及流量条件有效时滞。**Q2.2 于 2026-11-28 正式封存，不再开展新一轮模型搜索或实验。**

- 四阶段剂量关系显示 CO 流量由约 3.46 升至 43.18 m³/s 后趋于平台；CO₂ 流量在第二阶段约 18.34 m³/s 达到较高水平，随后下降。
- 单位耗氧边际排放随阶段变化：H0 下 CO 约 0.835、5.956、8.362、8.734；CO₂ 约 2.594、3.713、2.769、2.437 m³/m³。H1 仅作第五炉时基敏感性分析。
- 30 s 非重叠块的氧气增量与平均烟气流量呈中等但炉次异质的关联，不支持跨炉统一线性比例。
- 补充时滞检验改变窗口起点、流量估计口径，并进行已知 20 s 固定时滞注入。180 s CO₂ 的负向信号只在部分口径下出现；40 项主比较无一通过 Holm 校正。合成检验还显示中心差分流量可能制造虚假的负向关系。
- 因此，“高流量对应更短有效时滞”保留为未证实的物理假说，不宣布成立，也不宣布不存在。H0 为主情景，H1 为敏感性；不训练预测模型。

### Q2.3 有效响应时滞（已并入 Q2.2 并随之封存）

- 90 s 窗口下 CO、CO₂ 的有效时滞中位数均约 14 s；180 s 窗口下分别约 32 s、36 s。
- H0 180 s CO₂ 高低流量时滞中位差为 −33.5 s，H1 为 −33.0 s，但分别只有 4/5 炉能够形成炉内对比；窗口平移后没有足够窗口形成有效对比。
- 该时滞是反应、输运、取样和分析仪动态共同形成的有效响应时差，不是纯管道运输时间。十炉是独立实验单位，炉内窗口不是独立炉次。

[最终研究报告](reports/final/Q2_2_FINAL_RELATIONSHIP_STUDY.md) · [补充稳健性检验](reports/final/Q2_2_LAG_HYPOTHESIS_ROBUSTNESS.md) · [复算入口](experiments/q2_lag_robustness/) · [紧凑结果](results/q2/lag_hypothesis_robustness/)。

### Q2.3 第三轮：同炉同耗氧进度匹配（封存后探索性补充）

- 在同一炉次内，以累计耗氧进度为匹配变量，对低流量与高流量窗口进行无放回一对一匹配；90 s 为主窗口，180 s 为敏感性窗口，同时检查 Pearson/Spearman、不同流量分位及进度 caliper。
- H0、90 s、Pearson、25/75 分位、0.15 caliper 下，CO 与 CO₂ 各形成 7 对、覆盖 5/10 炉；炉次中位时滞差（高流量减低流量）分别为 +4 s、+8 s，95% 炉次 Bootstrap 区间分别为 [−4, +6] s、[−26, +20] s；两种气体均仅 1/5 炉呈负差。
- 0.05 和 0.10 进度 caliper 均无匹配；0.15 下配对进度差中位数为 0.134，说明主分析只是放宽容差后的阶段近似匹配，共同支持仍然不足。
- 因此第三轮仍未获得“高流量对应更短时滞”的稳健证据，也不能反向断言高流量会延长时滞。它是封存后的探索性 addendum，不改变原有冻结决定，不作因果解释。

[第三轮报告](reports/final/Q2_3_STAGE_MATCHED_LAG_ROUND3_ADDENDUM.md) · [复算入口](experiments/q2_stage_matched_lag/v1/) · [紧凑结果与 QA](results/q2/stage_matched_lag_2026-11-28/)。

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
| Q2.3 第三轮阶段匹配时滞 | [入口](experiments/q2_stage_matched_lag/v1/) | [结果](results/q2/stage_matched_lag_2026-11-28/) |
| Q2 时基审计 | [入口](experiments/q2_timebase_audit/) | [结果](results/q2/timebase_audit/) |
| Q3 B 炉 | [入口](experiments/q3_b_holdout/) | [结果](results/q3/b_holdout/) |
| Q3 梯度提升 | [入口](experiments/q3_boosted_curve/) | [结果](results/q3/boosted_curve/) |
| Q3 曲线迁移 | [入口](experiments/q3_curve_transfer/) | [结果](results/q3/curve_transfer/) |
| Q3 物理过程 | [说明](experiments/q3_physical_process/) | [结果](results/q3/physical_process/) |
| Q3 初始光谱气体预测 | [入口](experiments/q3_spectral_gas_prediction/) | [结果](results/q3/spectral_gas_prediction/) |
| Q3 V5 | [入口](experiments/q3_v5_causal_spectral/) | [结果](results/q3/v5_causal_spectral/) |

## 7. 下一步（建议优先级）

1. P0：将 Q2.1 与第二问最终关系研究整合进论文；明确 H0 主情景、H1/H2 局限及流量—时滞关系未获稳健支持的边界。
2. P0：统一 Q3 公平比较条件，决定每种气体的最终候选与退回模型，报告逐炉分数、误差和失败案例。
3. P0：撰写竞赛论文、摘要与控制页，引用冻结表和明确的探索性结果；不将报告目录中的历史版本当作最新交付。
4. P1：若可获得设备日志/时间戳，核验第五炉；否则保留不可识别结论及 H0/H1/H2 敏感性。不人为重建缺失原始观测。
5. P1：完成第二问图表、论文文字和全篇符号/单位的一致性检查；不再扩展第二问模型。

## 8. 后续提交约定

- 实验新增或变更时，记录问题编号、数据版本、信息集、验证协议、运行入口、输出路径、关键结果、证据限制和当前决策；新增专项同步更新第 6 节及首页导航。
- 已封存模型如需改变，新增版本并记录原因，不覆盖历史表格、不把探索性分数写成正式结论。全文报告优先于本总账的摘要数值。
- 结果大文件按既有 LFS/本地复算约定管理；历史资料留在 archive。默认批量脚本仅包含其显式注册的实验，不代表全部专项。
- 本次是文档和索引审查，未独立重新执行全部实验。
