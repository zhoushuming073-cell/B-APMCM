# 第二问时滞识别

复算有效滞后相关、负对照与 90 秒窗口留炉验证。光谱亮度诊断使用正式共同网格数据。

在根目录安装依赖后，运行 `python experiments/q2_delay_identification/run.py`。可传 `--output` 指定复算路径。

[已发表的实验结果](../../results/q2/delay_identification/)

迁移复核时，90 秒窗口的 190 行留炉分数与冻结结果一致。全炉配对块 bootstrap 的区间端点存在差异，汇总表的区间宽度中位数相差最多 1.99 s。历史发布表保存在 `tables/`，本次复算在 `rerun/`。
