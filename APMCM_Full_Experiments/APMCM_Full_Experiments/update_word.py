from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
import pandas as pd, os
src='/mnt/data/亚太杯.docx';out='/mnt/data/APMCM_Full_Experiments/亚太杯_完整数据实验更新.docx'
doc=Document(src)
# Update version line on cover if present
for p in doc.paragraphs:
    if '版本：V1.0' in p.text:
        p.text=p.text.replace('版本：V1.0','版本：V2.0')
    if '状态：建模方案已形成，仍待完整验证与定稿' in p.text:
        p.text='版本：V2.0｜2026年11月26日｜状态：十炉数据已完成第一、二问的主要探索实验；异常真值与时滞识别仍待验证'

def shade(cell,fill='D9EAF7'):
    tcPr=cell._tc.get_or_add_tcPr();shd=OxmlElement('w:shd');shd.set(qn('w:fill'),fill);tcPr.append(shd)
def set_cell(cell,text,bold=False):
    cell.text=str(text);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
    for p in cell.paragraphs:
        for r in p.runs:
            r.font.name='等线';r._element.rPr.rFonts.set(qn('w:eastAsia'),'等线');r.font.size=Pt(8.5);r.bold=bold

def add_table(headers,rows,widths=None):
    tab=doc.add_table(rows=1,cols=len(headers));tab.alignment=WD_TABLE_ALIGNMENT.CENTER
    for i,h in enumerate(headers):set_cell(tab.rows[0].cells[i],h,True);shade(tab.rows[0].cells[i])
    for row in rows:
        cells=tab.add_row().cells
        for i,v in enumerate(row):set_cell(cells[i],v)
    doc.add_paragraph('')
    return tab

def heading(text,level=1):
    p=doc.add_paragraph();p.style=doc.styles['Heading'+str(level)] if level>0 else doc.styles['Normal'];p.add_run(text)
    for r in p.runs:
        r.font.name='等线';r._element.rPr.rFonts.set(qn('w:eastAsia'),'等线')
    return p

def para(text,boldlead=None):
    p=doc.add_paragraph()
    p.paragraph_format.space_after=Pt(5)
    r=p.add_run(text);r.font.name='等线';r._element.rPr.rFonts.set(qn('w:eastAsia'),'等线');r.font.size=Pt(10)
    return p

# New appendix / version log
doc.add_page_break()
heading('15　完整十炉数据审计（V2.0新增）',1)
para('本轮数据来自用户上传的 Data.zip。十个 XLSB 工作簿均成功读取，每个工作簿一个工作表；所有序号连续，无光谱或过程变量缺失值。总计13,589个采样时点、27,816,683个光谱观测值。采样间隔按题面为0.5 s。原始工作簿未被修改。')
para('重要修正：前四炉6B系列共享完全相同的波长网格（340.535—1027.963 nm）；后六炉6C系列共享另一套网格（339.390—1023.303 nm）。两组不能按列号直接合并。跨组公共波长范围为340.535—1023.109 nm，共2030点；6C→公共网格→6C的线性插值往返MAE约28.6—40.2原始光强单位。组内模型应优先保留原生网格。')
para('全数据共有14,805个65535读数。每个饱和值在波长方向均至少有一个相邻饱和值；68.08%的饱和值前后两个时刻均未饱和，31.92%只有一侧时刻未饱和。共有6,735次累计耗氧相邻持平，没有负增量。CO₂零值146个、CO零值833个，均保留为观测值，不作为缺失值。')
audit=pd.read_csv('/mnt/data/APMCM_Full_Experiments/data_audit.csv')
add_table(['炉次','时点数','时长(s)','饱和点','耗氧持平步','CO₂零值','CO零值'],[[r.heat,int(r.rows),r.duration_s,int(r.saturated_65535),int(r.oxygen_flat_steps),int(r.CO2_zero_rows),int(r.CO_zero_rows)] for r in audit.itertuples()])

heading('16　第一问：全量数据实验',1)
heading('16.1　多尺度异常筛查',2)
para('每炉前55%时间段用于稳健尺度校准，后段作为测试，中间留40帧隔离。每炉注入144组人工异常事件，十炉共1440组、11,040个受干扰光强点；包括正负方向、5%/15%/40%幅度，以及单点、连续波长和连续时间异常。')
add_table(['检测器','平均点级召回率','平均事件级召回率'],[['时间Hampel','2.12%','5.42%'],['局部软联合残差','24.30%','58.06%'],['多尺度局部残差','88.24%','79.86%']])
para('多尺度模型按异常幅度的事件召回率为：5%约51.7%，15%约88.3%，40%约99.6%。它适合作为候选发现器，不适合作为自动删除规则。固定阈值在晚段的候选率平均约0.42%，高于校准段0.1%的预算，说明异常评分分布存在时间漂移。')
heading('16.2　双向模板复核与修复门控',2)
para('对筛查候选点，分别使用前、后相邻光谱，在局部峰肩上估计亮度校正。只有两个模板均可用、修正方向一致、最小偏差超过0.12个对数光强单位且偏差大于模板分歧的2倍时，才批准修复。全局门控参数在校准段人工注入实验中选定，再用于后段测试。')
para('十炉测试段总MAE宏平均由3.065降至2.217原始光强单位；人工注入点MAE由1713降至694；点级修复召回率47.8%，事件级修复召回率61.6%；总MAE在7/10炉改善。仍有未注入点被修改，且没有人工标注真值，因此该结果仅支持“保守修复情景”，不能视为已获得真实无异常光谱。')
para('将门控应用于全数据后，共筛出51,236个候选点，其中5,074点通过双向复核。所有原始值和掩码均保留。')
heading('16.3　跨炉共享光谱基底与峰值遮挡',2)
para('留一炉验证中，每次将一炉完全留出，用其余同标定组炉次训练基底；在测试炉中选择算法识别的局部峰候选，遮挡3、7、11个连续波长通道，仅用可见通道估计低秩系数。候选峰不是人工确认的物理发射线。')
add_table(['遮挡宽度','全炉混合秩12 MAE','同标定组秩12 MAE','同组原生网格秩18 MAE'],[['3通道','191.1','163.1','157.9'],['7通道','175.9','150.4','145.9'],['11通道','165.9','145.0','139.4']])
para('同标定组秩12基底在3/7/11通道测试中分别有8/10、9/10、9/10炉优于全炉混合基底。原生网格秩18在本轮遮挡实验中误差略低，但秩选择尚未完成独立嵌套验证；秩12保留为较保守基准。核心结论是两种标定网格不应强行共享同一个光谱基底。')
heading('16.4　饱和信号模拟与全炉修复情景',2)
para('实际14,805个饱和值中，14,373点获得一侧或双侧邻帧模板估计，432点暂未修复。另在十炉中各抽取20个未饱和峰候选，模拟50%和75%峰高截顶。邻帧模板的峰高相对误差中位数约19.6%，优于不修复和本轮简单删失高斯峰拟合；但峰面积误差并非在所有截顶程度下都更优。故饱和修复只作为条件估计，并保留下界与未修复标记。')
para('全炉暂定场景共修改19,447个光谱点：5,074个异常门控点和14,373个饱和模板点。它们保存在独立NPZ文件中，不覆盖原始Excel。十炉三维光谱图为显示而下采样，实际计算仍使用原始网格。')

heading('17　第二问：全量数据实验',1)
heading('17.1　累计耗氧与氧枪流量',2)
para('氧枪流量以累计耗氧曲线的导数表示。非负积分B样条以非负基函数系数表示潜在流量，并拟合其积分得到的累计耗氧。结点间距和平滑惩罚在内层时间块验证中选择，外层时间块独立评分。')
add_table(['模型','外层累计耗氧MAE(m³)','外层RMSE(m³)'],[['线性插值','2.361','2.859'],['PCHIP','5.278','6.137'],['非负积分样条','2.141','2.603']])
para('非负积分样条在7/10炉优于线性插值，但优势有限。该验证评价的是累计耗氧插值，不是瞬时流量真值，因为数据没有独立流量传感器标签。样条流量应与线性累计插值和更新区间平均斜率并列报告。')
para('累计耗氧共有6,735次相邻持平；正增量中位数约4.72—4.75 m³。九炉更新间隔中位数为1.0 s，6C02057为0.5 s。这与采样保持或计量刷新机制相符，但不能据此断言具体仪表机制。')
heading('17.2　CO₂、CO连续信号与累计排放',2)
add_table(['气体','线性MAE','PCHIP MAE','平滑样条MAE'],[['CO₂','0.0490','0.0504','0.1034'],['CO','0.1441','0.1869','0.2917']])
para('相同外层时间块验证下，线性插值是两种气体更稳健的简单基线；PCHIP作为形状保持对照。全炉累计排放按记录区间的梯形积分计算，不能外推记录窗口以外的排放。')
heading('17.3　氧量—烟气关系与过程阶段描述',2)
para('十炉Spearman相关系数中位数：累计耗氧与CO为0.817、与CO₂+CO为0.783、与CO₂为−0.155。差分后相应中位数约为0.210、0.216、−0.036。累计水平相关性受共同时间趋势影响，不能直接解释为因果关系。')
para('按时间等分为三段（不是经验证的冶金阶段），跨炉中位CO₂流量由17.06降至15.75再至13.07 m³/s；CO由3.82升至40.88再至43.12 m³/s；CO₂份额由0.791降至0.268再至0.239。该结果仅作描述性观察。')
heading('17.4　时滞识别与跨炉预测检验',2)
para('对每炉分别计算0—40 s内氧枪流量变化与烟气变化的互相关，并以循环移位作负对照，20个炉次×气体检验经BH校正后有6项q<0.05，但最佳时滞分布在0—40 s，最大绝对相关的中位数约0.123，部分峰值落在搜索边界。按窗口计算的流量与时滞Spearman相关约为CO₂ −0.023、CO −0.013，不支持稳定的“流量越大、时滞越短”单调关系。')
para('嵌套留一炉预测中，预测时点统一使用目标前5 s的信息。CO₂仅用历史CO₂/CO持续性基准MAE为0.311，加入氧枪当前或历史流量后反而变差（约0.344—0.352）。CO基准MAE为0.761，加入过程进度与历史气体后为0.734，加入当前氧枪流量为0.736；加入多段历史流量或滞后网格没有改善。该结果只表示预测增量，不是物理运输时滞的因果证明。')
para('九炉氧枪流量P10—P90重叠区约3.89—5.21 m³/s；6C02057明显偏高（约5.36/9.07/13.03 m³/s）。高流量支持主要由单炉提供，流量效应与炉次身份混杂。当前不应拟合并宣称已识别通用τ(q)函数。')

heading('18　当前推荐模型与未解决事项',1)
para('第一问：保留原生波长网格；多尺度时谱残差用于候选筛查，双向模板用于保守复核；共享基底按6B/6C标定组分别训练；秩12作为保守基准、秩18作为敏感性结果；饱和点保留下界和修复标记。')
para('第二问：非负积分样条作为氧枪流量主候选，线性累计插值和更新区间平均斜率作为对照；CO₂/CO采用线性插值主基线、PCHIP敏感性分析；时滞只报告分炉探索性结果，当前证据不足以支持单一τ(q)。')
para('仍需解决：①人工确认异常峰标签；②确认65535是否为仪器饱和上限及响应曲线；③将标定组效应与炉况效应分离；④获取独立瞬时氧枪流量或设备刷新日志；⑤提高高流量工况的跨炉覆盖；⑥对秩和门控进行严格嵌套验证；⑦第三问实时光谱—烟气预测尚未在本轮执行。')
para('完整复现代码、CSV结果、每炉连续曲线、三维图、候选清单和独立修复情景存放于 APMCM_Full_Experiments.zip；原始Data.zip未包含在结果包中，需由研究者另行保存。')

# Set metadata
doc.core_properties.title='亚太杯数学建模——完整十炉数据实验更新'
doc.core_properties.subject='2026 APMCM Problem B Questions 1 and 2 full-data experiments'
doc.core_properties.comments='V2.0 append-only research update; raw workbooks were not modified.'
doc.save(out)
print(out,os.path.getsize(out),'paras',len(doc.paragraphs),'tables',len(doc.tables))
