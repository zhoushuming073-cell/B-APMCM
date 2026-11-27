# 第三问曲线迁移、残差树与神经网络

在仓库根目录安装项目依赖并执行 Git LFS 拉取后运行。

    python experiments/q3_curve_transfer/run.py

复算写入 results/q3/curve_transfer/rerun。需要官方 Data.zip、十个 XLSB、process.csv 和 raw_*.npz。入口核对官方原文件哈希并核对缓存的过程通道。本次冻结实验另曾逐一核对十炉光谱与原始工作簿。

七种模型、三个情景、两种气体，共 420 个逐炉评分。模型超参数固定，主情景提前量在训练炉内选择。每个外层测试炉完整留出。独立炉次数为十。

[研究报告](../../reports/working/第三问异质性曲线机器学习与方差.md) · [完整预测](../../results/q3/curve_transfer/predictions.csv.gz) · [方差汇总](../../results/q3/curve_transfer/variance.csv)
