# Final Submission | 最终赛事交付

本目录为 2026 APMCM B 题最终交付的准备区，与研究记录、实验开发及历史归档分离。

## 赛氪上传清单

以实际提交页面为准。截止时间：2026-11-30 09:00（北京时间）。

| 文件 | 命名规则 | 格式 | 大小上限 | 状态 |
| --- | --- | --- | --- | --- |
| 英文论文（含 Summary Sheet） | B{TeamNumber}.pdf | PDF | 50 MB | 必交，待完成 |
| 承诺书（打印后全员手写签名） | B{TeamNumber}cns.pdf | PDF（亦支持图片） | 20 MB | 必交，待完成 |
| AI 使用说明 | B{TeamNumber}ai.pdf | PDF | 20 MB | 待完成 |
| 支撑材料 | B{TeamNumber}fj.zip | ZIP（亦支持 RAR/7z/PDF） | 200 MB | 选交，Q1/Q2 已归档，Q3 待封存 |

TeamNumber 为完整参赛队号（含 apmcm 前缀）。真实队号、签名及个人身份信息不应写入公开仓库。

## 支撑材料建议结构

```text
Supporting_Materials/
├── README.md
├── Mathematical_Models/
├── Source_Code/
├── Processed_Data/
├── Results/
├── Figures/
├── configs/
├── pyproject.toml
└── requirements.lock.txt
```

- 仅纳入与正式论文对应的模型、代码、必要处理数据、结果与复现说明。
- 不纳入赛题原始数据、Git LFS 指针、历史实验归档或临时文件。
- 确保 Q1–Q3 最终模型、图表、参数和运行入口与论文一致。
- AI 使用说明应如实记录实际工具、使用环节及人工核验；具体要求以赛事正式规范为准。
- 提交前检查论文匿名性、文件命名、PDF 可读性、文件大小及上传确认状态。
- 比赛期间不得将实际论文、签名承诺书或其他可能违反参赛纪律的材料公开；正式文件应在本地或获授权的私有环境准备。

## 状态

[Supporting_Materials](Supporting_Materials/) 已归档 Q1、Q2 冻结模型、源码、紧凑结果、最终验证及复现包（已排除直接原始数据副本）。Q3 待正式封存。最终图表、完整复现验收、论文、承诺书、AI 说明和最终 ZIP 仍需完成。既有研究资料见 [reports](../reports/) 与 [RESEARCH_STATUS.md](../RESEARCH_STATUS.md)。
