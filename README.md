# 2026 APMCM B 题 · 工作区导航

最后整理于 2026-11-26。查研究结论先读[实验进度](project/reports/实验进度.md)和[数学模型与实验证据](project/reports/2026%20APMCM%20B题——第三版数学模型与实验证据.md)。第三问开始建模前，先读[第三问项目指引](project/THIRD_QUESTION_PROJECT_GUIDE.md)。查第二问时滞结果优先看[90 秒窗口最终结果表](专题实验/Q2_Delay_Identification_Experiment/tables/final_delay_evidence_summary.csv)。

| 目录 | 用途与入口 |
|---|---|
| [2026 APMCM Problems](2026%20APMCM%20Problems/) | 赛事题面、论文规范与[更正版 Data.zip](2026%20APMCM%20Problems/2026%20APMCM%20Problem%20B/Data/Data.zip)。 |
| [project](project/) | 工作代码、测试、分版本数据、复核结果与持续更新的 Markdown 报告。[第三问项目指引](project/THIRD_QUESTION_PROJECT_GUIDE.md) · [当前数据审计](project/data/corrected_2026_11_26/audit.csv) |
| [APMCM_Full_Experiments](APMCM_Full_Experiments/APMCM_Full_Experiments/) | 全量实验原包，主工程从此读取十炉数据节点。[光谱曲面 PDF](APMCM_Full_Experiments/APMCM_Full_Experiments/figures/all_heats_spectral_surface.pdf) · [过程曲线 PDF](APMCM_Full_Experiments/APMCM_Full_Experiments/figures/all_heats_process_curves.pdf) |
| [专题实验](专题实验/) | 按研究问题保存输入、脚本和结果。第五炉时基、累计耗氧对照、时滞识别三个实验各自保留完整目录。 |
| [交付与图表](交付与图表/) | 已成稿的 V3 研究档案及第二问图表合集。 |

## 专题与交付入口

| 要查的问题 | 建议打开 |
|---|---|
| 第五炉为何保留两种采样间隔？ | [时间尺度审计](专题实验/第五炉时间尺度问题/report/APMCM_Timebase_Audit_Report.md) |
| 累计耗氧与烟气在 0 s、48 s 对齐下有何差异？ | [对照实验口径](专题实验/Q2_Aligned_Comparison/Q2_Aligned_Comparison/README_CN.txt) · [按炉留出汇总](专题实验/Q2_Aligned_Comparison/Q2_Aligned_Comparison/LOHO_summary.csv) |
| 供氧流量能否解释每个窗口的有效时差？ | [延迟识别报告](专题实验/Q2_Delay_Identification_Experiment/README_Q2_Delay_Experiment.md) · [逐炉误差](专题实验/Q2_Delay_Identification_Experiment/tables/90s_window_loho_models.csv) |
| 在哪里找可直接查看的报告和图片？ | [V3 研究档案 PDF](交付与图表/APMCM_B_V3_更新交付/APMCM_B_V3_Aligned_Research_Package/report/APMCM_B_V3_Aligned_Research_Record.pdf) · [可编辑 Word 版](交付与图表/APMCM_B_V3_更新交付/APMCM_B_V3_Aligned_Research_Package/report/APMCM_B_V3_Aligned_Research_Record.docx) · [第二问图表目录](交付与图表/Q2_Figures_Package/README_Q2_Figures.md) |

复现主工程时在 [project](project/) 下执行[进度表末尾的命令](project/reports/实验进度.md)。复现时滞专项时使用[包内运行说明](专题实验/Q2_Delay_Identification_Experiment/RUN_INSTRUCTIONS.md)。其中亮度诊断需要将完整实验包中的 [harmonized_spectra.npz](APMCM_Full_Experiments/APMCM_Full_Experiments/harmonized_spectra.npz) 放入该专项的 `inputs/`。

数据入口已复核。更正版 Data.zip 的 SHA-256 与 V3 记录一致；两个第二问专题及图表包各自附带的十份过程 CSV 与完整实验包逐文件哈希一致。

第三问正式结果待补。时滞专项的留炉误差单位为秒，五秒辅助预测允许使用历史烟气，这两套表均不能用来填报实时光谱预测精度。


## 从 GitHub 下载完整项目

原始附件 `Data.zip`、Excel 工作簿和 NPZ 数据由 Git LFS 管理。安装 Git LFS 后执行：

    git lfs install
    git clone https://github.com/zhoushuming073-cell/B-APMCM.git
    cd B-APMCM
    git lfs pull

更正版 `Data.zip` 的 SHA-256 为 `cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727`。
