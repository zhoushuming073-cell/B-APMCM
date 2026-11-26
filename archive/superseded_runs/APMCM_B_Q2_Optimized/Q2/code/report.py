"""Build the offline HTML research report from the recomputed tables."""
from html import escape
from reconstruct import ROOT,TAB,pd,np
from pathlib import Path

def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+str(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'

def fig(name,caption,height):
    return f'<figure><img style="max-height:{height}mm" src="../figures/{name}.png"><figcaption>{caption}</figcaption></figure>'

def make():
    pages=[]
    def page(kicker,title,body):
        n=len(pages)+1
        pages.append(f'<section><header>APMCM 2026　/　转炉炼钢过程<span>问题 B · 第二问</span></header><div class="eyebrow">{kicker}</div><h1>{title}</h1>{body}<footer>周树铭团队工作稿　·　第二问实验记录<span>{n:02d}　/　07</span></footer></section>')
    read=lambda n:pd.read_csv(TAB/(n+'.csv'))
    rec=read('reconstruction_holdout').query("scenario=='H0'").groupby(['channel','model'])[['mae','rmse','r2']].mean()
    scores=read('delay_folds').groupby(['scenario','window','gas','model']).mae.mean()
    rows=[]
    for k,lab in [('O','累计耗氧 / m³'),('CO2','CO₂ 流量 / m³/s'),('CO','CO 流量 / m³/s')]:
        a=rec.loc[(k,'selected')]
        rows.append([lab,f'{a.mae:.3f}',f'{rec.loc[(k,"linear baseline"),"mae"]:.3f}',f'{a.rmse:.3f}',f'{a.r2:.6f}'])
    page('01　研究结论','烟气曲线可以重构，流量效应仍难辨认',
         '<p class="lead">十炉、13,589 个原始时点。连续遮去短段后，CO₂ 与 CO 流量仍能较好恢复。把“供氧越快、延迟越短”写进模型后，两种气体在留出炉次上的平均误差均略升。</p>'
         +'<h2>观测缺口恢复精度</h2>'+table(['目标','本模型 MAE','线性插值 MAE','本模型 RMSE','本模型 R²'],rows)
         +'<p class="note">H0 主情景。每炉遮 12 段、共 144 个记录，每段 4 或 8 秒。先在每炉汇总被遮点，再对十炉等权平均。R² 描述恢复原始读数的精度，仪表自身误差和真实瞬时氧流量没有独立标定。</p>'
         +'<h2>跨炉结果改变了什么判断</h2><p>若流量确实提供了可迁移的延迟信息，加入流量参数应有助于解释未参与训练的烟气读数。固定时滞与流量非递增时滞的 MAE，CO₂ 分别为 <b>2.068 和 2.081</b>，CO 分别为 <b>6.847 和 6.865 m³/s</b>。这些差值很小，也没有跨设定稳定的改善。</p>'
         +'<div class="call"><b>第二问建议写法</b><br>给出十炉六类连续曲线。按流量报告局部时差的中位数与四分位区间。保留满足题面方向的输运模型作条件情景，标明哪些秒数来自模型约束、哪些来自观测对齐。</div>'
         +'<h2>原始数据复核</h2><p>附件 1 的 CO₂ 与 CO 最高分别为 21.74 和 59.04 m³/s，均有变化。十炉同时纳入。三个过程通道从 XLSB 独立提取，与前序实验的三个输入通道逐元素核对，最大绝对差为零。附件表头称 concentration，积分时采用题面第 3 页明确给定的体积流量单位。</p>')
    page('02　连续重构','从累计表读数恢复连续氧流量',
         '<h2>耗氧表保留端点与单调性</h2><p>用清理后的原始读数连成单调折线 L(t)。在首末点作中心对称延拓，取宽度 w 的滑动平均。这样首末累计量不变，流量连续且非负，分段积分能回到账面净增量。</p>'
         +'<div class="eq">Ô(t) = (1/w) ∫<sub>t−w/2</sub><sup>t+w/2</sup> L(u) du<br>q̂(t) = [L(t+w/2) − L(t−w/2)] / w</div>'
         +'<p>训练炉比较 4、8、16 秒平均窗及三个宽度的保形样条。H0 的十个外层折中，9 折选择 4 秒平均窗，1 折选择 8 秒。另保留原始折线和原始节点 PCHIP 对照。折线导数在拐点跳变，仅作为累计量误差基线。</p>'
         +'<h2>烟气节点与连续排放量</h2><p>分别对 CO₂、CO 做孤点筛查，再将邻近记录按时间窗口汇成节点。节点时刻取观测时间的均值，节点值取流量均值，保留首末观测，用保形三次 Hermite 插值连接 [1]。两种气体在 H0 的十折均选中 2 秒节点。非负节点间保形，避免负排放。</p>'
         +'<div class="eq">Ĝ<sub>g</sub>(t) = PCHIP(t<sub>j</sub>, Ḡ<sub>g,j</sub>)(t)　　 Ê<sub>g</sub>(t) = ∫<sub>0</sub><sup>t</sup> Ĝ<sub>g</sub>(u) du</div>'
         +'<p>对分段多项式解析积分。与原始读数的梯形积分相比，十炉末端累计排放最大相对差为 CO₂ 0.0163%、CO 0.0207%。这是重构前后总量的一致性检查。</p>'
         +'<h2>异常规则能检出什么</h2><p>耗氧表如出现回退，比较回退两侧样点相对邻点插值的偏差，优先剔除偏差较大者。本批原始数据无回退。气体要求疑点在相邻两点之间形成孤立尖峰或凹陷、两侧相互接近，且偏差超过 max(0.2, 6 × 1.4826 × 局部 MAD)。原始数据共标记 3 条 CO、0 条 CO₂；疑点和处理值可逐行查看。</p>'
         +'<p>人工注入每档每通道 447 个孤立扰动。±4 m³/s 的 CO₂、CO 分别检出 438、407 个；±1 m³/s 时降至 357、219 个。耗氧 ±15 与 ±60 m³ 均为 447/447。此检验不覆盖连续漂移或持续跳变，也不能判断这 3 个原始候选是否为设备故障。</p>')
    a=read('dose_response')
    dt=a.pivot_table(index='oxygen_net',columns='outcome',values='value',aggfunc='median')
    rows=[[f'{i:,.0f}',*[f'{dt.loc[i,k]:,.2f}' for k in ['CO2','CO','E_CO2','E_CO']]] for i in [500,1500,2500,3000]]
    page('03　供氧与烟气','同一耗氧进度下，CO 与 CO₂ 走势不同',
         '<p>将每炉起点扣除，令 X(t) = Ô(t) − Ô(0)。横轴统一为已投入的氧量。按 X 查回原曲线，在十炉间取中位数。气体仍处于传感器时钟上，因此这里描述同一记录时刻的关系。</p>'
         +fig('dose_response','蓝线为各炉，橙线为同一耗氧量处可用炉次的中位数。下表四个位置均覆盖十炉。',92)
         +table(['净耗氧 / m³','CO₂ / m³/s','CO / m³/s','累计 CO₂ / m³','累计 CO / m³'],rows)
         +'<h2>累计曲线接近完全相关，局部关系弱得多</h2><p>逐炉 Spearman 相关再取中位数，净耗氧对 CO₂、CO 瞬时流量分别为 −0.155、0.817，对两种累计排放量均约为 1。累计量共同随时间增长，会机械地推高相关。</p><p>将每炉切成不重叠的完整 30 秒段，比较段内新增氧量与气体平均流量，CO₂、CO 的相关中位数降至 0.289、0.335；与两种新增排放量的相关分别为 0.295、0.335。后两组近似来自同一积分恒等关系。不能据此识别供氧的因果弹性。</p>')
    page('04　时滞建模','用烟气原始读数检验对齐方案',
         '<p>同一炉内，氧的变化要经过反应、混合、输运和仪器响应才能体现在烟气读数上。其他钢厂工艺研究也分别记录了取样管输运与分析仪响应 [3]。本题没有示踪脉冲、仪器阶跃响应及独立流量真值，可检验的是信号之间的有效时差。</p>'
         +'<h2>带流量的时间坐标</h2><div class="eq">z(q) = clip((q − 5)/5, −1, 1)，　v(q) = exp[β z(q)]<br>C(t) = ∫<sub>0</sub><sup>t</sup> v(q(u)) du，　s(t) = C<sup>−1</sup>(C(t) − d)，　τ(t) = t − s(t)</div>'
         +'<p>t 是气体读数时刻，s(t) 是与它配对的供氧时刻，d 是 5 m³/s 参考流量下的秒数。v 恒正保证反解唯一。β = 0 退化为固定时滞；β &gt; 0 在恒流局部给出 τ(q) = d / v(q)。流量变化时计算积分再反解。反方向 β &lt; 0 也作为对照。</p>'
         +'<h2>控制吹炼阶段，再比较供氧对齐</h2><p>传感器时间统一取 2 秒网格。阶段向量含截距、t/600、净耗氧/3000 及其在 0.15、0.40、0.70、0.90 的正部折线项。将 q(s)/5 及其平方加入阶段向量，用岭回归拟合原始烟气读数。炉等权、截距免罚、惩罚系数 0.001。纯阶段模型不含流量项，零时滞模型使用同时刻流量。</p>'
         +'<p>时滞试验单独固定 8 秒中心差商作为输入，4、16 秒检查敏感性。此处两端将差商截在有效时间域内；该处理与第二页端点反射的绘图曲线略有不同。它是离线对齐，中心差商利用邻近未来耗氧，不用于声称在线预测。</p>'
         +'<h2>整炉嵌套留出</h2><p>外层留一炉，其余九炉再逐炉留一选择参数 [2]。网格为 d = 0, 6, …, 60 秒，β 取 0、±ln(1.5)、±ln(2)。测试时同一时刻、同一原始气体读数、同一阶段基线。按炉先算 MAE，再等权平均。v 最小为 0.5，故回看最多 120 秒；主比较统一从第 120 秒开始。吹炼最初 120 秒不在该评分的适用范围内。</p>'
         +'<div class="call">另用十炉真实供氧轨迹人工生成烟气响应，在已知响应函数、加噪条件下，固定与递减两类 24 秒参考延迟各重复 40 次，均找回植入网格点。它检查反演和搜索实现，在真实炉况中仍需留炉检验。</div>')
    rows=[]
    for scen in ['H0','H1','H0_without_5']:
        for gas in ['CO2','CO']:
            rows.append([{'H0':'H0','H1':'H1','H0_without_5':'剔除第 5 炉'}[scen],gas.replace('2','₂'),*[f'{scores.loc[(scen,8,gas,m)]:.3f}' for m in ['stage','zero','fixed','decreasing','reverse']]])
    page('05　检验结果','强制时滞递减，留出炉次没有稳定获益',
         table(['情景','气体','仅阶段','零时滞','固定','非递增','反向允许'],rows)
         +'<p class="note">原始烟气流量 MAE，单位 m³/s，越低越好。固定、非递增和反向允许各自在内层选参，均允许 d = 0。剔除第五炉后重新执行全部内外层拟合。</p>'
         +fig('delay_paired','每根柱是同一留出炉的“非递增 − 固定”MAE。负值表示流量依赖模型改善。',68)
         +'<p>H0 下，CO₂ 的炉均差为 +0.0126 m³/s，2 炉改善、4 炉持平、4 炉变差；CO 为 +0.0179，4 炉改善、3 炉持平、3 炉变差。将已有的十个配对误差按炉重采样，95% 描述区间分别为 [−0.0002, 0.0320]、[−0.0291, 0.0758]。训练集在各外层折间重叠，此区间不代表独立新厂区的覆盖率。</p>'
         +'<h2>改动一个设定，最优时差就可能改变</h2><p>H0 将氧流量差商改为 4 或 16 秒，CO 的“递减 − 固定”误差差值分别为 −0.098 和 −0.001 m³/s；CO₂ 对应 +0.008、+0.011。H1 将第五炉时长从 387 改为 774 秒，平均供氧从 9.139 降至 4.569 m³/s，其余炉次时基不变。第五炉重构 CO 总排放由 11,030.6 变为 22,063.4 m³，绝对时间尺度仍缺独立记录佐证。</p>')
    b=read('empirical_flow_bins').query("scenario=='H0'")
    rows=[]
    for band,sub in b.groupby('band',sort=False):
        q=sub.iloc[0]
        cells=[]
        for gas in ['CO2','CO']:
            r=sub[sub.gas==gas].iloc[0]
            cells.append(f'{r.lag_median:.0f} [{r.lag_q25:g}, {r.lag_q75:g}]')
        rows.append([{'low':'低','middle':'中','high':'高'}[band],f'{q.q_min:.2f}–{q.q_max:.2f}',int(q.n_windows),*cells])
    page('06　分流量报告','同一流量区间内，候选延迟相差十余秒',
         '<h2>数据直接给出的局部对齐结果</h2><p>每炉取不重叠的 90 秒窗口，使用平滑后的氧流量和气体变化率。每个候选滞后都用相同的前 30 个输入点、步长 2 秒，在 0–30 秒中找绝对相关最大处。将全部窗口依其前 60 秒平均氧流量分成三组。窗口包含多个反应，峰位仅是候选时差。</p>'
         +table(['流量组','实际覆盖 / m³/s','窗口数','CO₂ 时差 / s','CO 时差 / s'],rows)
         +'<p class="note">单元格为中位数 [第 25, 第 75 百分位]，不是置信区间。最高组内流量分布偏斜，中位氧流量仅 5.47 m³/s。0 秒与 30 秒的边界峰未剔除。</p>'
         +fig('local_flow_bins','主时基与第五炉替代时基下，中流量组均出现较短的局部峰位。',58)
         +'<p>在炉内同时对流量与峰位扣除线性时间趋势，残差回归斜率为 CO₂ −0.50、CO −0.45 秒/(m³/s)。按整炉抽样 5,000 次的区间为 [−4.23, 3.24] 与 [−3.71, 2.86]。各炉对斜率的贡献取决于其残差变异和窗口数。</p>'
         +'<h2>题面给定递减方向时的条件计算</h2><p>用十炉留一的均值选取最终展示参数：H0 的 CO₂ 为 d = 6 秒、β = ln 2；CO 为 d = 60 秒、β = ln 1.5。取恒流 4、5、6 m³/s，CO₂ 的条件时滞为 <b>6.89、6.00、5.22 秒</b>，CO 为 <b>65.07、60.00、55.33 秒</b>。CO 的参考延迟顶到搜索上限，不能定为实际输运时间。H1 重新标定后，两气体的参考延迟分别变为 60、12 秒。</p>'
         +'<p class="note">局部互相关搜索的是 90 秒窗口内的短期峰；输运模型针对跨炉气体水平误差选参，60 秒上限施加在参考 d 上。二者估计目标不同，不应将两组秒数平均。</p>')
    lag=read('global_lags').query("scenario=='H0'");em=read('emissions').query("scenario=='H0'").set_index('heat')
    rows=[]
    for h,g in lag.groupby('heat'):
        e=em.loc[h]
        cells=[]
        for gas in ['CO2','CO']:
            v=g[g.gas==gas].iloc[0]
            cells.append(f'{v.lag_s:.0f} / {v.bh_q:.3f}')
        rows.append([h,f'{e.duration:g}',f'{e.net_O:,.1f}',f'{e.E_CO2:,.1f}',f'{e.E_CO:,.1f}',*cells])
    page('07　炉次核查','十炉排放总量与全炉有效时差',
         table(['炉号','时长 / s','净耗氧 / m³','CO₂ 总量 / m³','CO 总量 / m³','CO₂ 秒 / q','CO 秒 / q'],rows)
         +'<p class="note">末两列为绝对互相关峰位 / BH 校正值。相同的 2 秒网格，8 秒氧差商；对两序列用标准差 3 秒的高斯平滑后取梯度，扫描 0–60 秒。比较时固定输入段；用远离零位移至少 90 秒的循环错位作负对照，每个错位重新搜索全部滞后。每情景校正 20 个检验。非平稳信号的循环错位仅是近似对照 [4]。</p>'
         +'<p>CO₂ 在 7/10 炉通过 q ≤ 0.10，CO 在 0/10 炉通过。全炉峰位中位数为 35 与 27 秒，中位峰强 |r| 为 0.216 与 0.139。CO₂ 有部分可重复的时间对齐迹象；对 CO，表中秒数只宜保留为探索性结果。通过错位对照也不能证明峰位精确。</p>'
         +'<h2>参考资料</h2><ol><li>Fritsch &amp; Butland (1984). A Method for Constructing Local Monotone Piecewise Cubic Interpolants. <i>SIAM J. Sci. Stat. Comput.</i> DOI 10.1137/0905021.</li><li>Roberts et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. <i>Ecography.</i> DOI 10.1111/ecog.02881.</li><li>Bathla (2011). Carbon Prediction Model Using the Offgas Analysis at ArcelorMittal Indiana Harbor No. 2 Steel Producing BOF Shop. <i>Iron &amp; Steel Technology.</i> 完整题名与原文入口见文献目录。</li><li>Pennsylvania State University. <i>STAT 510, Lesson 9.1.</i> Pre-whitening as an Aid to Interpreting the CCF.</li></ol>'
         +'<p class="note">数据来源为 2026 APMCM Problem B 及附件 1–10。十张分炉六联图、原始过程列、候选参数全表、逐折误差与预测、异常点及合成检验明细均随包附上。</p>')
    css="""@page{size:A4;margin:0}*{box-sizing:border-box}body{margin:0;background:#ecefeb;color:#1c303a;font:10.4pt/1.55 'Noto Sans CJK SC','Microsoft YaHei',sans-serif}section{width:210mm;height:297mm;position:relative;padding:13mm 16mm 15mm;background:white;margin:0 auto 5mm;break-after:page}section:last-child{break-after:auto}header,footer{font-size:8pt;letter-spacing:.3px;color:#58717c;border-bottom:1px solid #cedbdd;padding-bottom:3mm}header span,footer span{float:right}footer{position:absolute;bottom:8mm;left:16mm;right:16mm;border:0;border-top:1px solid #cedbdd;padding:1.7mm 0 0}.eyebrow{font-weight:700;color:#008083;margin:8mm 0 1.5mm;font-size:9pt}h1{font-size:21pt;line-height:1.3;margin:0 0 5mm;letter-spacing:-.8px}h2{font-size:12.6pt;margin:4mm 0 1.5mm;color:#156266}p{margin:2mm 0}.lead{font-size:12pt;color:#234851}.note,figcaption{color:#48616c;font-size:8.9pt;line-height:1.5}.call,.eq{background:#eaf3f2;border-left:3px solid #14817e;padding:3mm 4mm;margin:3mm 0}.eq{font:12pt/1.8 'DejaVu Serif','Noto Sans CJK SC',sans-serif}table{border-collapse:collapse;width:100%;margin:2mm 0;font-size:9.3pt}td,th{padding:1.3mm 1.3mm;text-align:right;border-bottom:1px solid #d9e2e5;font-variant-numeric:tabular-nums}td:first-child,th:first-child{text-align:left}th{background:#e9f1f2;color:#21515a}figure{margin:3mm 0}img{display:block;margin:auto;max-width:100%;object-fit:contain}figcaption{text-align:center;margin-top:1mm}li{margin:1mm 0}ol{padding-left:6mm;font-size:9.5pt} @media print{body{background:white}section{margin:0}*{-webkit-print-color-adjust:exact;print-color-adjust:exact}}"""
    f=ROOT/'report'/'Q2_report.html'
    f.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>第二问｜转炉供氧与烟气</title><style>'+css+'</style><body>'+''.join(pages)+'</body></html>',encoding='utf-8')

if __name__=='__main__':
    make()
