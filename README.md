# 2026 APMCM B 题

转炉火焰光谱连续重构、供氧与烟气时滞分析。正式数据使用更正版的十个炉次，共 13,589 帧。第三问已完成十炉、双时基、19 个模型族的探索性嵌套留炉实验；光谱增量尚未稳定证实。

| 目录 | 内容 |
| --- | --- |
| [configs](configs/) | 官方数据位置、第一问重构及各第二问实验的输入和复算目录 |
| [docs](docs/) | 赛题、官方规范和方法与验证协议 |
| [data](data/) | 原始附件、提取数据、处理后数据和审计 |
| [src/apmcm](src/apmcm/) | 数据读取、谱面重构和烟气实验代码 |
| [experiments](experiments/) | 各专项的运行入口与输入说明 |
| [results](results/) | 第一、二、三问的冻结实验结果 |
| [reports](reports/) | 研究文档、已成稿报告和提交材料分区 |
| [archive](archive/) | 迁移前的原始交付包和历史工作目录 |

## 数据与复现

先安装 Git LFS。以下命令在仓库根目录运行。依赖锁按 Windows 上的 Python 3.13 解析。

```powershell
git lfs pull
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python scripts/verify_data.py
python -m pytest
```

```powershell
python scripts/run_all.py --help
python scripts/run_all.py q2_optimized
python scripts/run_all.py q3_spectral_gas_prediction
python scripts/run_all.py
python scripts/build_figures.py
```

批量复算结果位于各问题的 `rerun/`，图形重绘位于 `results/rebuilt_figures/`。运行入口使用本目录中的配置和正式数据。校验失败时可使用 `git lfs pull` 取得实际二进制文件。原始附件 SHA-256 与解压后的逐炉文件均可核对。

- 正式数据与审计 [Data.zip](data/raw/official/Data.zip) · [SHA-256](data/metadata/checksums.sha256) · [审计表](data/metadata/audit.csv)
- 第一问 [光谱曲面汇编](results/q1/figures/all_heats_spectral_surface.pdf) · [窄峰遮挡指标](results/q1/tables/q1_peak_mask_fit_metrics.csv)
- 第二问 [优化实验报告](reports/final/q2_optimized/Q2_report.pdf) · [90 秒窗口结果](results/q2/delay_identification/tables/final_delay_evidence_summary.csv) · [时基审计](reports/final/timebase_audit/APMCM_Timebase_Audit_Report.pdf)
- 第三问 [模型及实测报告](reports/working/Q3_SPECTRAL_GAS_MODEL_REPORT.md) · [实验入口](experiments/q3_spectral_gas_prediction/) · [冻结结果](results/q3/spectral_gas_prediction/)
- 研究进度 [实验进度](reports/working/实验进度.md) · [数学模型与实验证据](reports/working/2026%20APMCM%20B题——第三版数学模型与实验证据.md) · [第三问项目指南](reports/working/THIRD_QUESTION_PROJECT_GUIDE.md) · [第三问研究记录](reports/working/第三问研究.docx)
