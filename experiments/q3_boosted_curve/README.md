# 第三问梯度提升追加实验

在仓库根目录安装依赖并执行 Git LFS 拉取后运行。

    python experiments/q3_boosted_curve/run.py
    python experiments/q3_boosted_curve/plot.py --source results/q3/boosted_curve/rerun

实验复用 curve_transfer 的数据、配准和评分逻辑，以固定复杂度直方图梯度提升修正共享曲线。运行还会读取上一轮冻结的 folds.csv，校验相同提前量和目标点数并保存比较表。模型参数固定，提前量用训练炉内部验证选择。输出默认写入 results/q3/boosted_curve/rerun。

[详细报告](../../reports/working/第三问梯度提升追加实验20个R2.md) · [20 个 R²](../../results/q3/boosted_curve/twenty_r2.csv) · [20 张图](../../results/q3/boosted_curve/figures/README.md)
