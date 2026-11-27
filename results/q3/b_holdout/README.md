# B 炉单独留出：冻结结果

- selected_summary.csv：三轮由训练炉内层选定的模型；V2/V3 外部指标保留六位小数，详见原始表。
- v2/、v3/：原始全部候选 B 指标（108/308 行）、九炉逐折指标（972/2772 行）、两气体内层排名（各54/154行），已与本地实验表 Git blob 哈希核对一致。
- V1 原始全表、分类元数据、逐点预测及预测图保存在原交付 Q3_B_Trial.zip；本目录保留其选定模型汇总。
- V2/V3 逐点预测、完整缓存及运行日志保存在 Q3_B_Optimization_V2_V3.zip，未将大体积预测或中间数据重复放入 Git。按实验入口可复算。

原交付 ZIP SHA-256：
Q3_B_Trial.zip: 16a3f5f1b52cd9ffb3ade63cfb109155e4dc2725114d94ca02f15f5ff5877ed6
Q3_B_Optimization_V2_V3.zip: cf0c9b433b49a4916260abbcdc355529d965f2aa6472ef3041bb2936b7561ee0

这些是 B-only 探索结果，不是十炉正式泛化性能。详见[报告](../../../reports/working/Q3_B_HOLDOUT_MODEL_REPORT.md)。
