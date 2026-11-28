# Q3 物理化学与测量语义审计 V2

本专项优先核实 CO/CO₂通道究竟是浓度还是体积流量，并用氧平衡检验流量解释，再量化末段气体组成变化。

## 数据与协议

- 正式 Data.zip，SHA-256：cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727
- 十炉、13,589 行，H0：0.5 s/行。
- 5.5 s 时间轴中值滤波；氧枪供氧速度由累计氧差分后再中值滤波估计。
- 条件流量假设检验：若 CO/CO₂为同基准 m³/s，计算 ∫(CO₂ + 0.5 CO)dt / 累计氧枪供氧。
- 尾段比较终点前 60 s 的起点 10 s 与终点 10 s 中位数；第五、第七炉按既有终端断崖规则截断。

## 关键结果

- 条件流量假设下，十炉 COx 氧当量/氧枪累计供氧比为 2.93–6.53，中位数 5.77；该解释与氧平衡不相容。
- 尾段端点比较：CO₂通道上升、CO 通道下降、COx 通道和下降、表观 CO₂份额上升，均为 10/10 炉。
- 表观 CO₂份额不是完整 PCR；COx 通道和不是绝对碳排放量。
- 第五炉末段供氧速度估计增加 82.4%，其余九炉首尾窗口变化在 ±5% 内。

## 文件

- measurement_audit.py：从正式 Data.zip 读取十炉原始 XLSB，输出审计表。
- measurement_semantics_audit.csv：逐炉条件氧平衡比及通道统计。
- tail_composition_oxygen_audit.csv：逐炉尾段组成与供氧速度端点比较。
- run_manifest.json：数据哈希、样本量和关键汇总。

## 复算

~~~bash
PYTHONPATH=src python -m apmcm.q3.measurement_audit --data-zip data/raw/official/Data.zip --output-dir results/q3/physical_chemistry_v2
~~~

## 证据边界

题面“concentration”与“m³/s”描述冲突，原始 XLSB 表头也有中英文及 YQXH 变体。体积分数是合理解释，但需工厂数据字典确认。没有总烟气流量、H₂、H₂O、O₂、N₂和浴内 C/FeO/温度等数据，不能识别绝对脱碳速率、完整 PCR 或各反应速率。十炉多轮探索不构成独立外部验证。