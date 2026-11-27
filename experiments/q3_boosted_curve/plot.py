"""Plot all 20 held-out heat/gas predictions from the saved boosting run."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[2]


def main(source):
    points = pd.read_csv(source/"predictions.csv.gz").query("model=='boost_fusion' and scenario=='H0_q2'")
    scores = pd.read_csv(source/"folds.csv").query("model=='boost_fusion' and scenario=='H0_q2'")
    heats = sorted(points.heat.unique())
    assert len(heats)==10 and len(scores)==20
    font_names = {f.name for f in font_manager.fontManager.ttflist}
    for name in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC"]:
        if name in font_names:
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 11
    out = source/"figures"
    out.mkdir(exist_ok=True,parents=True)
    records = []
    limits = {g: float(np.ceil(points.query("gas==@g")[["observed","forecast"]].max().max()/5)*5)
              for g in ["CO2","CO"]}

    def draw(ax, heat, gas, compact=False):
        data=points.query("heat==@heat and gas==@gas").sort_values("detector")
        s=scores.query("heat==@heat and gas==@gas").iloc[0]
        assert len(data)==s.n
        r2=1-np.sum((data.observed-data.forecast)**2)/np.sum((data.observed-data.observed.mean())**2)
        assert abs(r2-s.r2)<1e-7
        label="CO2" if gas=="CO2" else "CO"
        ax.plot(data.detector,data.observed,color="#173d68",lw=1.55,label="实际值",zorder=2)
        ax.plot(data.detector,data.forecast,color="#d76224",lw=1.5,ls="--",label="预测值",zorder=3)
        ax.set_title(f"附件 {heats.index(heat)+1} · {heat} · {label}     R² = {s.r2:.4f}",fontsize=10 if compact else 13,loc="left",pad=12)
        ax.set_ylim(0,limits[gas])
        ax.set_xlim(data.detector.min(),data.detector.max())
        ax.set_xlabel("检测时刻（s）")
        ax.set_ylabel("瞬时流量（m³/s）")
        ax.grid(alpha=.18)
        ax.spines[["top","right"]].set_visible(False)
        ax.legend(loc="upper right",frameon=False,ncol=2,fontsize=8 if compact else 10)
        if not compact:
            ax.text(0,-.24,f"整炉留出验证  |  提前 {s.lead:.0f} 秒预测  |  MAE = {s.mae:.3f} m³/s  |  RMSE = {s.rmse:.3f} m³/s",
                    transform=ax.transAxes,fontsize=10,color="#495565")
        return s

    with PdfPages(out/"全部20组预测与实际.pdf") as pdf:
        for i,h in enumerate(heats,1):
            for gas in ["CO2","CO"]:
                fig,ax=plt.subplots(figsize=(10,4.6))
                fig.subplots_adjust(left=.085,right=.98,top=.86,bottom=.25)
                s=draw(ax,h,gas)
                file=out/f"{i:02d}_{h}_{gas}.png"
                fig.savefig(file,dpi=170,facecolor="white")
                pdf.savefig(fig,facecolor="white")
                plt.close(fig)
                records.append(dict(attachment=i,heat=h,gas=gas,r2=s.r2,lead_s=s.lead,image=file.name))
    for gas in ["CO2","CO"]:
        fig,axes=plt.subplots(5,2,figsize=(15,17),layout="constrained")
        for h,ax in zip(heats,axes.flat):
            draw(ax,h,gas,True)
        fig.suptitle(f"{'CO2' if gas=='CO2' else 'CO'} · 十炉留出预测与实际",fontsize=18)
        fig.savefig(out/f"{gas}_overview.png",dpi=150,facecolor="white")
        plt.close(fig)
    pd.DataFrame(records).to_csv(out/"index.csv",index=False)
    print(f"Saved {len(records)} individual PNGs, two overview PNGs and a 20-page PDF to {out}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=ROOT/"results/q3/boosted_curve")
    main(parser.parse_args().source)
