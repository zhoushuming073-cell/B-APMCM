# 第三问：实时光谱—烟气软测量

以更正版原始十炉数据直接执行嵌套留一炉验证；H0 为主，H1 为第五炉时基敏感性情景。所有模型仅使用发布时刻及以前的光谱和累计耗氧，无在线烟气输入。

仓库根目录：

    git lfs pull
    python -m pip install -r requirements.lock.txt
    python -m pip install --no-deps -e .
    python experiments/q3_spectral_gas_prediction/run.py

可用 --archive 指定已核对 SHA-256 的 ZIP，--output 指定输出路径。默认复算写入 results/q3/spectral_gas_prediction/rerun/，不覆盖冻结表。输出完整逐点预测、逐炉指标、超参数、时滞、配对对照、负对照及四张十炉预测图。

[模型报告](../../reports/working/Q3_SPECTRAL_GAS_MODEL_REPORT.md) · [冻结结果](../../results/q3/spectral_gas_prediction/) · [模型代码](../../src/apmcm/q3/model.py)

本轮为两阶段探索性验证，未证明光谱稳定增量或通用物理运输时间，不应把候选集成模型描述为已完成现场验收。
