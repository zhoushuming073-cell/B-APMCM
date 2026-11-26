# 2026 APMCM B · 第二问优化实验

**先读 report/Q2_report.pdf**（7 页）。report/Q2_report.html 可离线打开，保留 figures 与 report 的相对位置。每个炉次各一张六联图，文件名为 figures/heat_附件号_炉号.png。

## 在指定 Windows 工作目录运行

将 ZIP 解压至 `C:\Users\周树铭\Desktop\APMCM codexB`，得到其中的 `Q2` 目录。已附第二问完整的三个原始过程通道，共 13,589 行。复现全部数值实验无需寻找或再次下载 2,047 列光谱。

已安装 Python 3.12 且可联网安装依赖时，双击 `Q2\run_windows.cmd`。脚本在项目目录建立独立虚拟环境，读取 data/process.csv，重算 tables 和 figures。也可从 PowerShell 执行：

```powershell
Set-Location -LiteralPath 'C:\Users\周树铭\Desktop\APMCM codexB\Q2'
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
$env:OPENBLAS_NUM_THREADS = '1'
.\.venv\Scripts\python.exe code\run.py
.\.venv\Scripts\python.exe code\validate.py
```

如要从原始工作簿重新提取，使用 `.\.venv\Scripts\python.exe code\extract.py "..\Data.zip"`（替换为原包实际位置），之后再运行 run.py。extract.py 针对本题 XLSB 单表、指定列与数值类型；可逐炉核查 data/sources.csv 的 SHA-256。完整光谱和原始 Data.zip 不重复打包。报告给出的曲线对齐、缺口恢复均是离线实验。Windows 命令按标准 Python 启动器编写，实验与校验在 Linux、Python 3.12 中实际运行。

## 如何读结果

| 文件 | 内容 |
|---|---|
| data/process.csv、sources.csv | 附件号、炉号、行号、原始累计耗氧 O、原始 CO₂ 与 CO 体积流量；各 XLSB 的校验值与原表头 |
| tables/curves.csv.gz | 两种时基的逐点六条曲线，O、q、CO2、CO、E_CO2、E_CO。另附原始列。pandas.read_csv 可直接打开 gzip |
| tables/reconstruction_holdout.csv | 留炉选参后，在被遮原始读数上的 MAE、RMSE、R²；每炉每通道共 144 个评估点 |
| tables/reconstruction_cv.csv、reconstruction_choice.csv | 各候选、掩码折、选择结果；linear 是原始折线，clean_linear 是筛除疑点的折线，pchip 为未聚合的保形插值，aW 为均宽 W 秒的滑动平均，kW 为 W 秒分窗的保形插值 |
| tables/curve_audit.csv、spike_flags.csv、spike_positive_control.csv | 全观测贴合度、孤点清单、已知人工尖峰检出试验 |
| tables/emissions.csv、associations.csv、dose_response.csv | 逐炉累计排放、分段关联、不同净耗氧量对应的两种流量与两种累计量 |
| tables/local_windows.csv、empirical_flow_bins.csv、flow_slope.csv | 局部候选峰和流量分箱、控制炉内时间趋势后的斜率与按炉重采样区间 |
| tables/global_lags.csv、global_profiles.csv | 各炉两种气体的全局对齐、多重检验校正、每个候选滞后的带符号相关 |
| tables/delay_folds.csv、delay_predictions.csv.gz | 嵌套留炉误差、参数和主带宽下逐点预测。所有对照使用同一测试时刻 |
| tables/delay_profiles.csv、delay_calibration.csv、delay_by_flow.csv | 固定候选参数的留炉风险、全样本展示参数、按指定流量计算的条件延迟与支持炉数 |
| tables/paired_delay.csv、delay_positive_control.csv、validation_mass.csv | 配对差值、时滞反演阳性对照、每炉氧量积分一致性 |
| figures/ | 十炉六联图以及流量分箱、相关剖面、逐炉误差、第五炉时基和耗氧关联图 |
| references/ | 参考文献、DOI/出处、每篇具体用途、BibTeX |

H0 = 全部按题面每行 0.5 秒；H1 = 仅第五炉每行 1 秒；H0_without_5 = 剔除第五炉，整套时滞预测重新校准。H1 不是数据订正。d 是 5 m³/s 参考流量下的延迟。经验峰位和给定单调规律后的条件延迟是两种不同的量，详见报告第 6 页。

## 复现口径

实验入口 code/run.py 调用重构、时间对齐、合成检验、制图；code/validate.py 检查端点、非负与单调性、积分守恒、十炉六联图及模型间测试行一致性。code/report.py 从表格生成 HTML 版式，叙述及关键数字对应此次冻结结果；可在浏览器打印为 PDF。固定随机种子 20261126。数值依赖见 requirements.txt。data/archive_sha256.txt 对应原始 Data.zip。ARTIFACT_SHA256.txt 用于核验解压后的交付文件。
