"""Held-out heat prediction and variance accounting with no gas feedback."""
import argparse
from pathlib import Path
import platform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from threadpoolctl import threadpool_limits

from apmcm.common.paths import PROC, PROCESS_INPUT, RESULTS
from apmcm.data.verify import verify, digest
from apmcm.q3.curve_transfer import bands, prepare, Experiment, GASES


def load_inputs():
    verify()
    process = pd.read_csv(PROCESS_INPUT)
    raw = {}
    provenance = []
    for h, frame in process.groupby("heat"):
        file = PROC/f"raw_{h}.npz"
        with np.load(file, allow_pickle=False) as a:
            f = frame.sort_values("row")
            np.testing.assert_array_equal(f.row, a["sequence"])
            np.testing.assert_allclose(f[["O","CO2","CO"]],np.column_stack([a["oxygen"],a["gas"]]),rtol=0,atol=1e-8)
            x = a["x"]
            if x.shape != (len(f),2047) or not np.isfinite(x).all() or not np.all(np.diff(a["wl"])>0):
                raise ValueError(f"Spectra failed audit for {h}")
            raw[h] = (a["oxygen"],a["gas"],bands(a["wl"],x))
        provenance.append(dict(heat=h,rows=len(f),attachment=int(f.attachment.iloc[0]),sha256=digest(file)))
    if len(raw)!=10 or sum(len(v[0]) for v in raw.values())!=13589:
        raise ValueError("Invalid heat coverage")
    return raw,pd.DataFrame(provenance)


def summarize(rows):
    keys = ["scenario","model","gas"]
    folds,summary = [],[]
    for (scenario,model,gas,h), f in rows.groupby(keys+["heat"]):
        y, yhat = f.observed.to_numpy(), f.forecast.to_numpy()
        err = yhat-y
        folds.append(dict(scenario=scenario,model=model,gas=gas,heat=h,lead=f.lead.iloc[0],n=len(f),
                          r2=1-err@err/np.sum((y-y.mean())**2),mae=np.mean(abs(err)),
                          rmse=np.sqrt(np.mean(err**2)),bias=err.mean(),residual_var=err.var(ddof=1)))
    folds=pd.DataFrame(folds)
    for (scenario,model,gas),f in folds.groupby(keys):
        summary.append(dict(scenario=scenario,model=model,gas=gas,macro_r2=f.r2.mean(),
                            r2_var=f.r2.var(ddof=1),r2_sd=f.r2.std(ddof=1),min_r2=f.r2.min(),
                            macro_mae=f.mae.mean(),macro_rmse=f.rmse.mean(),
                            mean_within_residual_var=f.residual_var.mean(),
                            between_bias_var=f.bias.var(ddof=1),heats_above_09=int((f.r2>.9).sum())))
    return folds,pd.DataFrame(summary)


def run(out):
    raw, manifest = load_inputs()
    out.mkdir(exist_ok=True,parents=True)
    manifest.to_csv(out/"inputs.csv",index=False)
    points = []
    for scenario,timebase,scheme in [("H0_q2","H0","q2"),("H0_30s","H0","30s"),("H1_30s","H1","30s")]:
        data={h:prepare(*v,1. if h=="6C02057" and timebase=="H1" else .5) for h,v in raw.items()}
        exp=Experiment(data)
        for h, a in data.items():
            train=tuple(i for i in data if i!=h)
            for g,gas in enumerate(GASES):
                lead=exp.lead(train,g,scheme)
                end,idx,pred=exp.predict(train,h,g,lead)
                for name,forecast in pred.items():
                    points.append(pd.DataFrame(dict(scenario=scenario,heat=h,model=name,gas=gas,
                                  lead=lead,issue=a.t[idx],detector=a.t[end],
                                  observed=a.target[end,g],forecast=forecast)))
        exp.clear()
        print(f"Completed {scenario}",flush=True)
    points=pd.concat(points,ignore_index=True)
    folds,summary=summarize(points)
    points.to_csv(out/"predictions.csv.gz",index=False,float_format="%.10g",compression={"method":"gzip","mtime":0})
    folds.to_csv(out/"folds.csv",index=False)
    summary.to_csv(out/"variance.csv",index=False)
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True,layout="constrained")
    for gas,ax in zip(GASES,axes):
        for model,group in folds.query("scenario=='H0_q2' and gas==@gas").groupby("model"):
            ax.plot(group.heat,group.r2,marker=".",label=model)
        ax.set_ylabel(gas+" held-out R²")
        ax.axhline(0,color="gray",lw=.7)
    axes[0].legend(ncol=3,fontsize=8)
    axes[1].tick_params(axis="x",rotation=30)
    fig.savefig(out/"fold_r2.png",dpi=170)
    plt.close(fig)
    print("Python",platform.python_version(),"sklearn",sklearn.__version__,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=RESULTS/"q3/curve_transfer/rerun")
    with threadpool_limits(1):
        run(p.parse_args().output)
