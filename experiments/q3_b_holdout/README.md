# Q3 B炉专项探索（V1–V3）

固定附件2（6B03206）为外部测试炉，A、C–J 为训练炉；H0 0.5秒/行。每版本以九炉内部留炉平均 R² 选候选，B 不参与拟合。B 曾被用于探索及反复评估，不能当作全新确认集。结果不覆盖原十炉实验。

运行（项目根目录，已安装依赖）：

    git lfs pull
    python experiments/q3_b_holdout/run.py --versions 1 2 3

原实验使用的精确中间特征缓存可从此前交付的 Q3_B_Trial.zip 或 Q3_B_Optimization_V2_V3.zip 解压取得；若保留原缓存，推荐：

    python experiments/q3_b_holdout/run.py --cache /path/to/q3_eda --versions 1 2 3

未指定 cache 时入口从正式 ZIP 重新构造物理谱带均值、95%分位与中位对比度。该新仓库适配入口未完成端到端复跑；原始缓存版模型源码仅调整导入与路径，已核对内容。原始源码和中间特征缓存见交付包。复算默认写入 results/q3/b_holdout/rerun/v1–v3，不覆盖冻结表。

[研究报告](../../reports/working/Q3_B_HOLDOUT_MODEL_REPORT.md) · [冻结结果](../../results/q3/b_holdout/)
