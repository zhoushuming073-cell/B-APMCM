from __future__ import annotations
import json,sys
from pathlib import Path
import pandas as pd

def table(df,cols,rename=None,fmt=None):
    z=df[cols].copy()
    if rename:z=z.rename(columns=rename)
    if fmt:
        for c,f in fmt.items():
            if c in z:z[c]=z[c].map(lambda x: f.format(x) if pd.notna(x) else '')
    return z.to_markdown(index=False)

def write_report(results_dir,report_path):
    r=Path(results_dir);read=lambda n:pd.read_csv(r/n)
    actual=read('actual_lag_summary.csv');bins=read('actual_flow_bins.csv');match=read('matching_summary_old_new.csv')
    reg=read('flow_lag_regression_old_new.csv');noiseless=read('noiseless_self_check_summary.csv')
    recovery=read('synthetic_recovery_equal_heat.csv');zero=read('zero_response_equal_heat.csv')
    fpr=read('zero_response_false_positive_summary.csv');scale=read('synthetic_scale_diagnostics_summary.csv')
    qa=json.loads((r/'qa_manifest.json').read_text());cfg=json.loads((r/'config.json').read_text())
    h0=actual[(actual.metric=='pearson')&(actual.scenario=='H0')&(actual.window_s==90)].copy()
    h0180=actual[(actual.metric=='pearson')&(actual.scenario=='H0')&(actual.window_s==180)].copy()
    h1=actual[(actual.metric=='pearson')&(actual.scenario=='H1')].copy()
    b0=bins[(bins.metric=='pearson')&(bins.scenario=='H0')&(bins.window_s==90)].copy()
    b0['flow_bin']=pd.Categorical(b0.flow_bin,['low','middle','high'],ordered=True);b0=b0.sort_values(['gas','flow_bin'])
    m15=match[(match.caliper==.15)].copy();m10=match[(match.caliper==.10)].copy()
    rec90=recovery[recovery.window_s==90].groupby(['gas','true_lag_s','snr'],as_index=False).agg(
        equal_heat_bias_s=('equal_heat_bias_s','mean'),equal_heat_mae_s=('equal_heat_mae_s','mean'),
        equal_heat_rmse_s=('equal_heat_rmse_s','mean'),within_2s=('equal_heat_within_2s_rate','mean'),
        boundary_rate=('equal_heat_boundary_rate','mean'),direction_accuracy=('equal_heat_direction_sign_accuracy','mean'))
    rec180=recovery[recovery.window_s==180]
    nl=noiseless.groupby(['gas','window_s','true_lag_s'],as_index=False).agg(
        median_estimate=('median_estimated_lag_s','mean'),within_2s=('within_2s_rate','mean'),
        boundary_rate=('boundary_rate','mean'),direction_accuracy=('direction_sign_accuracy','mean'))
    z90=zero[zero.window_s==90];z180=zero[zero.window_s==180]
    fprtab=fpr.copy()
    overlap=read('flow_progress_overlap_summary.csv')
    scale_by_gas=scale.groupby('gas',as_index=False).agg(actual_residual_sd=('actual_residual_sd','mean'),synthetic_residual_sd=('synthetic_residual_sd','mean'),
        actual_residual_acf1=('actual_residual_acf1','mean'),synthetic_residual_acf1=('synthetic_residual_acf1','mean'),
        injected_noise_acf1=('injected_noise_acf1','mean'),realized_snr=('realized_snr','mean'),
        actual_missing_fraction=('actual_missing_fraction','mean'),synthetic_missing_fraction=('synthetic_missing_fraction','mean'))
    def rows(df,columns,labels=None):return table(df,columns,labels)
    lines=[]
    lines += ['# Q2.3 最终验证与封存报告','',f"- 基线 GitHub commit：`4333a0bd9267cb899867b9e4bc0abc286267e8f0`",f"- 数据输入 SHA-256：`{qa['source_data_zip_sha256']}`",f"- 重建曲线 SHA-256：`{qa['curve_sha256']}`",f"- 随机种子：`{qa['seed']}`；每个随机条件重复：`{qa['simulation_replicates_per_condition']}` 次",'', '## 首页结论：三个问题的直接回答','',
      '### 1. 两种气体在不同流量下测到什么？','',
      'H0、90 秒主窗口的绝对互相关峰位中位数均为 14 秒。按全体窗口的氧枪流量三分位分组，CO 的低/中/高流量组中位峰位为 13/12/16 秒，CO₂ 为 13/14/14 秒。没有出现稳定的随流量增加而缩短的单调关系。峰位是表观时移，不是已确认的物理输运时间。','',
      rows(b0,['gas','flow_bin','cut_low_m3_s','cut_high_m3_s','flow_min_m3_s','flow_max_m3_s','n_windows','n_heats','lag_median_s','lag_q25_s','lag_q75_s'],{'gas':'气体','flow_bin':'流量组','cut_low_m3_s':'全局Q33','cut_high_m3_s':'全局Q67','flow_min_m3_s':'组内流量最小','flow_max_m3_s':'组内流量最大','n_windows':'窗口数','n_heats':'炉数','lag_median_s':'峰位中位数(s)','lag_q25_s':'Q25(s)','lag_q75_s':'Q75(s)'}), '',
      '### 2. 在什么信号条件下估计器可信？','',
      '无噪声、已知延迟自检能够验证时间方向、单位和索引：90 秒窗口下，0/10/20/30 秒注入延迟的中位估计均回到真值附近，±2 秒命中率约 80%–83%；30 秒真值恰好位于搜索上界，约 79% 的窗口命中上界。180 秒窗口把最大搜索范围同时扩展到 60 秒后，30 秒注入不再处于边界，理想条件下命中率接近 100%。','',
      '但在真实十炉供氧输入、真实阶段趋势和其他炉次经验残差块构成的半合成背景中，90 秒窗口的 MAE 约为 6–13 秒，±2 秒命中率通常只有约 28%–45%。无响应控制的 90 秒峰位中位数为 CO 14.4 秒、CO₂ 14.8 秒，与真实数据的 14 秒非常接近。故 14 秒不能作为可信的绝对延迟测量。','',
      '180 秒窗口的零响应控制更严重：原循环移位检验在 100 个零响应面板中，CO 有 37% 的面板至少出现一个 BH q≤0.1 的窗口，CO₂ 为 75%；窗口级 q≤0.1 比例分别约 15.2% 和 32.1%。因此 180 秒版本的显著窗口数不宜作为可靠证据。','',
      '### 3. 十炉数据对流量独立效应支持到什么程度？','',
      '同炉、相近耗氧进度的高低流量匹配，在容差 0.10 下为 0 对；放宽到 0.15 后各气体 7 对，来自相同 5 炉。修正进度代表点后，匹配差仍为 CO +4 秒、CO₂ +8 秒，炉次 Bootstrap 区间均跨零。线性/二次阶段调整下系数会改变符号，置信区间均跨零。数据没有支持独立负向流量效应，也不足以证明其不存在。','',
      '## 一、锁定口径与实际时滞结果','',
      '主分析：H0、90 秒窗口、2 秒网格、最大绝对 Pearson 相关；氧枪输入使用窗口前 2/3（60 秒），烟气使用完整 90 秒；正滞后表示氧枪变化在前、烟气检测在后。对氧枪输入和烟气分别线性去趋势，AR(1) 系数由氧枪序列估计并用于两者预白化。平局时取最小滞后。','',
      '官方题目说明浓度栏的 CO₂、CO 数值按瞬时体积流量处理，单位 m³/s；累计耗氧量单位 m³。H0 时间步长为 0.5 秒；H1 为既有替代时间轴敏感性情景（6C02057 使用 1 秒，其余 0.5 秒）。','',
      '### H0 主窗口及 H0 180 秒敏感性','',
      rows(pd.concat([h0.assign(时间窗='90s'),h0180.assign(时间窗='180s')]),['时间窗','gas','n_windows','n_heats','lag_median_s','lag_q25_s','lag_q75_s','positive_peak_n','negative_peak_n','zero_lag_n','upper_bound_n','bh_q_le_0_1_n','bh_denominator','min_empirical_p'],{'时间窗':'窗口','gas':'气体','n_windows':'窗口数','n_heats':'炉数','lag_median_s':'峰位中位数(s)','lag_q25_s':'Q25(s)','lag_q75_s':'Q75(s)','positive_peak_n':'正相关峰','negative_peak_n':'负相关峰','zero_lag_n':'零滞后峰','upper_bound_n':'上界命中','bh_q_le_0_1_n':'BH q≤0.1','bh_denominator':'分母','min_empirical_p':'最小经验p值'}), '',
      '90 秒窗口每个窗口只有 13 个可用对照移位，最小经验 p 值为 1/14≈0.0714。0/71 表示当前检验体系下没有窗口通过 BH 阈值，不等于证明不存在关联。180 秒窗口有 28 个可用移位，最小经验 p 值为 1/29≈0.0345。','',
      '### H1 替代时间轴敏感性','',rows(h1,['window_s','gas','n_windows','lag_median_s','lag_q25_s','lag_q75_s','bh_q_le_0_1_n','bh_denominator'],{'window_s':'窗口(s)','gas':'气体','n_windows':'窗口数','lag_median_s':'峰位中位数(s)','lag_q25_s':'Q25(s)','lag_q75_s':'Q75(s)','bh_q_le_0_1_n':'BH q≤0.1','bh_denominator':'分母'}),'',
      '## 二、已知延迟恢复实验','',
      '注入层：原始 Data.zip 提取出的十炉气体原始观测时间戳层，先用 Q2.1 冻结氧枪输入构造响应，再生成合成原始烟气观测，随后经过同一 k2 局部均值/PCHIP 重构、2 秒重采样、线性去趋势、AR(1) 预白化、窗口筛选与寻峰。没有把实测烟气曲线平移后当作真值。','',
      '合成响应：y_syn(t)=T_gas(t)+β·sign·Q_fluct(t−τ)+ε_block(t)。T_gas 为 60 秒居中移动平均；Q_fluct 为氧枪流量减去同一 60 秒移动平均；β 按目标信噪比设定，使注入响应与经验噪声的原始层标准差比为 1 或 0.5。ε_block 来自其他炉次同种气体去趋势残差的连续移动块，块长由十炉残差 ACF 首次降至 1/e 的中位滞后决定：CO 11.5 秒（23 点），CO₂ 12.25 秒（25 点）。经验残差可能包含真实工艺波动，故称经验背景而非纯测量噪声。','',
      '经验噪声块长由残差 ACF 首次低于 1/e 的滞后确定：CO 11.5 秒（23 点），CO₂ 12.25 秒（25 点）。下面比较实测与半合成序列的残差波动、自相关和缺测率；半合成残差的短期自相关略低于实测，故该背景只能视为近似经验背景。','',
      rows(scale_by_gas,['gas','actual_residual_sd','synthetic_residual_sd','actual_residual_acf1','synthetic_residual_acf1','injected_noise_acf1','realized_snr','actual_missing_fraction','synthetic_missing_fraction'],{'gas':'气体','actual_residual_sd':'实测残差SD','synthetic_residual_sd':'合成残差SD','actual_residual_acf1':'实测残差ACF(1)','synthetic_residual_acf1':'合成残差ACF(1)','injected_noise_acf1':'注入噪声ACF(1)','realized_snr':'实际注入SNR','actual_missing_fraction':'实测缺测率','synthetic_missing_fraction':'合成缺测率'}),'',
      '### 90 秒主场景（按方向和信噪比平均）','',
      rows(rec90,['gas','true_lag_s','snr','equal_heat_bias_s','equal_heat_mae_s','equal_heat_rmse_s','within_2s','boundary_rate','direction_accuracy'],{'gas':'气体','true_lag_s':'真值(s)','snr':'目标SNR','equal_heat_bias_s':'Bias(s)','equal_heat_mae_s':'MAE(s)','equal_heat_rmse_s':'RMSE(s)','within_2s':'±2s命中率','boundary_rate':'上界命中率','direction_accuracy':'方向符号识别率'}), '',
      '### 180 秒扩展搜索敏感性','',
      rows(rec180,['gas','true_lag_s','direction','snr','equal_heat_mae_s','equal_heat_within_2s_rate','equal_heat_boundary_rate','equal_heat_direction_sign_accuracy'],{'gas':'气体','true_lag_s':'真值(s)','direction':'方向','snr':'目标SNR','equal_heat_mae_s':'MAE(s)','equal_heat_within_2s_rate':'±2s命中率','equal_heat_boundary_rate':'上界命中率','equal_heat_direction_sign_accuracy':'方向符号识别率'}), '',
      '30 秒真值在 90 秒主窗口恰好等于最大搜索上界；扩展对照使用 180 秒窗口和 60 秒最大搜索滞后，因此窗口长度与搜索上界同时改变，不能把差异单独归因于搜索上界。','',
      '### 无噪声自检','',rows(nl,['gas','window_s','true_lag_s','median_estimate','within_2s','boundary_rate','direction_accuracy'],{'gas':'气体','window_s':'窗口(s)','true_lag_s':'真值(s)','median_estimate':'中位估计(s)','within_2s':'±2s命中率','boundary_rate':'上界命中率','direction_accuracy':'方向符号识别率'}),'',
      '## 三、零响应对照与显著性校准','',
      '零响应控制保留阶段趋势和经验背景残差，但响应系数严格设为零。每种气体、每个随机条件 100 次背景重采样；重复是模拟次数，不是增加独立炉次数。','',
      rows(zero,['gas','window_s','n_heats','equal_heat_median_peak_lag_s','equal_heat_peak_q25_s','equal_heat_peak_q75_s','equal_heat_bh_false_positive_rate'],{'gas':'气体','window_s':'窗口(s)','n_heats':'炉数','equal_heat_median_peak_lag_s':'零响应峰位中位数(s)','equal_heat_peak_q25_s':'Q25(s)','equal_heat_peak_q75_s':'Q75(s)','equal_heat_bh_false_positive_rate':'q≤0.1比例'}),'',
      rows(fpr,['gas','window_s','n_replicates','replicates_with_any_discovery','replicates_with_any_q_rejection_rate','mean_q_rejections_per_replicate','mean_window_q_rejection_rate'],{'gas':'气体','window_s':'窗口(s)','n_replicates':'重复数','replicates_with_any_discovery':'至少一次拒绝的面板数','replicates_with_any_q_rejection_rate':'面板任一q拒绝比例','mean_q_rejections_per_replicate':'每面板平均拒绝数','mean_window_q_rejection_rate':'窗口级q拒绝比例'}),'',
      '90 秒零响应控制的峰位分布与真实 14 秒中位数高度重叠；这说明当前估计器在该经验背景下可能产生相似峰位。不能据此分解真实 14 秒中“多少秒来自噪声”，也不能将经验背景解释为纯噪声。180 秒移位检验在零响应控制下出现较高拒绝比例，说明其离散 p 值与窗口依赖会导致校准问题，180 秒显著窗口数不宜作为确认性证据。','',
      '## 四、流量效应识别审计','',
      '### 旧/新进度代表点的匹配对照','',rows(m15,['progress_definition','gas','caliper','n_pairs','n_heats','median_progress_gap','max_progress_gap','median_flow_diff_m3_s','median_heat_lag_diff_s','bootstrap_heat_median_ci_low','bootstrap_heat_median_ci_high','n_heats_high_flow_shorter','n_heats_high_flow_longer'],{'progress_definition':'进度口径','gas':'气体','caliper':'容差','n_pairs':'匹配对数','n_heats':'炉数','median_progress_gap':'进度差中位数','max_progress_gap':'最大进度差','median_flow_diff_m3_s':'流量差中位数(m³/s)','median_heat_lag_diff_s':'按炉中位滞后差(s)','bootstrap_heat_median_ci_low':'炉次Bootstrap下限','bootstrap_heat_median_ci_high':'炉次Bootstrap上限','n_heats_high_flow_shorter':'高流量更短炉数','n_heats_high_flow_longer':'高流量更长炉数'}),'',
      '容差 0.10 下无匹配；容差 0.15 下两种气体均 7 对、5 炉。旧口径为整窗中点进度，新口径为氧枪输入前 2/3 区间的平均累计耗氧进度。修正后进度差中位数由约 0.134 降至约 0.119，但配对数和滞后差方向没有实质改变。两种气体共享相同工况窗口，不得把 7+7 对视为 14 对独立证据。','',
      '共同支持诊断：十炉中有 5 炉的低/高流量进度 min–max 区间相交，但实际低流量窗口没有落入这些交集区间（0 个），高流量窗口有 20 个。这说明区间范围相交并不等于拥有充分的实际观测重叠，0.15 容差匹配本质上仍是较宽松的近邻比较。','',
      rows(overlap,['progress_definition','n_heats','heats_with_progress_overlap','total_low_windows_inside_overlap','total_high_windows_inside_overlap','median_overlap_progress_width'],{'progress_definition':'进度口径','n_heats':'炉数','heats_with_progress_overlap':'进度范围相交炉数','total_low_windows_inside_overlap':'交集内低流量窗口','total_high_windows_inside_overlap':'交集内高流量窗口','median_overlap_progress_width':'交集宽度中位数'}),'',
      '### 连续阶段调整回归（第四轮原线性/二次规格）','',rows(reg,['progress_definition','gas','degree','boundary_mode','n_windows','n_heats','beta_s_per_within_heat_sd','cluster_se','ci_low','ci_high'],{'progress_definition':'进度口径','gas':'气体','degree':'阶段多项式阶数','boundary_mode':'峰位样本','n_windows':'窗口数','n_heats':'炉数','beta_s_per_within_heat_sd':'系数(s/炉内流量SD)','cluster_se':'炉次聚类SE','ci_low':'95%下限','ci_high':'95%上限'}),'',
      '线性与二次调整的流量系数会改变符号，置信区间均跨零。新旧进度代表点的系数差异很小，不改变结论。共同支持不足仍然是主要识别瓶颈。','',
      '## 五、第二轮复现与版本审计','',
      f"- 第二轮修复版复算了 40 项主检验（5 种氧枪流量构造 × 2 个时间轴 × 2 种气体 × 2 种窗口），Holm 校正后通过 0 项。\n- 仓库归档的 `key_actual_results.csv` 只有 24 条精选记录，并非完整 40 行。此次复算与归档 24/24 条记录中所有已归档字段一致（部分记录本身未保存全部检验字段）；其余 16 项原始完整结果此前未归档，无法进行历史逐项比较。\n- 90 秒最小可达经验 p=1/14；180 秒最小可达经验 p=1/29。\n- 旧 `prewhitened_delay` 版本在 180 秒下报告 CO 0/33、CO₂ 12/33；正式 `relationship_study_final/lag_summary.csv` 报告 CO 16/33、CO₂ 0/33。该翻转与预处理/代码版本同时变化，仓库缺少逐步骤单因素复验，故保留为未定位的历史版本差异；正文只引用正式 `lag_summary.csv`，旧表标为历史敏感性材料。",'',
      '## 六、输入与复现审计','',
      f"- 原始 `Data.zip` SHA-256 与仓库冻结清单一致：`{qa['source_data_zip_sha256']}`。\n- 提取的 `process.csv` Git blob SHA-1 与 main 中该文件一致：`{qa['process_csv_github_blob_sha1']}`；但冻结 manifest 中记录的 process CSV SHA-256 与当前文件内容不一致，已在差异审计表中保留，不静默覆盖历史记录。\n- 冻结 Q2.1 氧枪流量 profile 的 20 炉次摘要字段与 main 中正式 `flow_summary.csv` 全部一致；其 SHA-256 为 `11f3363e...`。\n- 重建曲线共 27,178 行，10 炉、H0/H1 两套时间轴，联合键唯一；重建曲线 SHA-256 与正式 `relationship_study_final` 输入指纹 `9bdc1fe5...` 完全一致。\n- 与第三轮输入指纹 `c89a8e34...` 的 27,178 行数据逐键匹配；数值最大绝对差仅约 `4.55e-13`，均值约 `5.34e-15`，属于浮点序列化量级差异。\n- 第四轮引用的 `f7280a7e...` 历史曲线文件未随仓库保存，故无法做该版本逐行比较。\n- H0 原始数据缺失值为 0；合成输入也没有引入缺测。",'',
      '### 冻结曲线输入核验','',
      rows(read('flow_profile_archive_comparison.csv'),['scenario','heat','numeric_fields_match','numeric_fields_compared','fit_success_match','all_fields_match'],{'scenario':'时间轴','heat':'炉次','numeric_fields_match':'匹配数值字段','numeric_fields_compared':'数值字段总数','fit_success_match':'拟合状态一致','all_fields_match':'全部一致'}),'',
      rows(read('curve_reconstruction_archive_comparison.csv'),['rows_current','rows_round3','unique_key_current','unique_key_round3','unmatched_rows','max_abs_numeric_difference','mean_abs_numeric_difference','current_matches_relationship_study_hash','round3_reference_hash_matches_round3_manifest'],{'rows_current':'当前行数','rows_round3':'第三轮行数','unique_key_current':'当前联合键唯一','unique_key_round3':'第三轮联合键唯一','unmatched_rows':'未匹配行','max_abs_numeric_difference':'最大绝对数值差','mean_abs_numeric_difference':'平均绝对数值差','current_matches_relationship_study_hash':'正式关系研究哈希一致','round3_reference_hash_matches_round3_manifest':'第三轮哈希一致'}),'',      '100 次半合成恢复实验使用的曲线文件为第三轮输入指纹 `c89a8e34...`；当前从官方输入与冻结 Q2.1 profile 重建的正式曲线为 `9bdc1fe5...`。两者 27,178 行逐键匹配，最大数值差为 4.55e-13。使用固定随机种子做的一次完整模拟窗口对照，2,536/2,536 行的延迟、峰值和符号结果完全相同；真实数据 570 个窗口的所有关键字段也完全相同。故本轮 100 次模拟结果在机器精度差异下可视为数值等价，但报告保留两个输入指纹，不声称 100 次均已在 9bdc 字节版本上重跑。','',
      (rows(read('actual_window_equivalence_check.csv'),['field','rows_compared','all_equal','max_abs_difference'],{'field':'字段','rows_compared':'窗口数','all_equal':'完全一致','max_abs_difference':'最大差值'}) if (r/'actual_window_equivalence_check.csv').exists() else '当前单次复跑目录未附历史窗口级对照表；正式滞后汇总仍与归档正式表逐项核对。'),' ',
      (rows(read('simulation_curve_equivalence_check.csv'),['field','rows_c89','rows_9bdc','all_equal','max_abs_difference'],{'field':'字段','rows_c89':'第三轮窗口数','rows_9bdc':'正式窗口数','all_equal':'完全一致','max_abs_difference':'最大差值'}) if (r/'simulation_curve_equivalence_check.csv').exists() else '当前独立复跑直接使用正式 9bdc 曲线，未使用第三轮 c89 曲线做模拟等价性对照。'),' ',
      '## 七、最终判断与停止规则','',
      '**最终状态：带明确限制封存。**','',
      '可以交付的结论是：十炉数据下 CO、CO₂ 的表观互相关峰位随窗口尺度变化；90 秒主窗口的经验中位峰位为 14 秒，但零响应经验背景也能产生相近峰位，故不将其解释为物理输运时间。按低、中、高流量分组未观察到稳定递减；同炉、近进度匹配和阶段调整均未支持独立负向流量效应，且共同支持不足。','',
      '不允许写成：高流量必然缩短时滞；CO/CO₂ 已确认存在统一 14 秒物理延迟；第四轮证明流量效应不存在；半合成恢复成功等于独立工业输运测量。','',
      '本轮完成后停止 Q2.3 模型搜索。只保留代码、数据版本、文稿和图表的工程修正，不再新增候选模型。','']
    Path(report_path).write_text('\n'.join(lines),encoding='utf-8')

if __name__=='__main__':
    write_report(sys.argv[1],sys.argv[2])
