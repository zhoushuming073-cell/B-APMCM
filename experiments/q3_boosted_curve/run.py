"""Run a fixed-complexity gradient-boosting benchmark and save every heat score."""
import argparse
from pathlib import Path
import sys
import platform
import json
import numpy as np
import pandas as pd
import sklearn
from threadpoolctl import threadpool_limits
from apmcm.common.paths import ROOT, RESULTS
from apmcm.q3.curve_transfer import prepare, Experiment, GASES
from apmcm.q3.boosted_curve import predict

sys.path.insert(0, str(ROOT))
from experiments.q3_curve_transfer.run import load_inputs, summarize


def run(out):
    raw, manifest = load_inputs()
    data = {h:prepare(*v,.5) for h,v in raw.items()}
    experiment = Experiment(data)
    out.mkdir(exist_ok=True, parents=True)
    manifest.to_csv(out/"inputs.csv",index=False)
    points = []
    for h,a in data.items():
        train=tuple(k for k in data if k!=h)
        for g,gas in enumerate(GASES):
            lead=experiment.lead(train,g,"q2")
            end,idx,forecasts=predict(experiment,train,h,g,lead)
            for name,forecast in forecasts.items():
                points.append(pd.DataFrame(dict(scenario="H0_q2",heat=h,gas=gas,model=name,
                              lead=lead,issue=a.t[idx],detector=a.t[end],
                              observed=a.target[end,g],forecast=forecast)))
        print(h,flush=True)
    points=pd.concat(points,ignore_index=True)
    folds,summary=summarize(points)
    points.to_csv(out/"predictions.csv.gz",index=False,float_format="%.10g",compression={"method":"gzip","mtime":0})
    folds.to_csv(out/"folds.csv",index=False)
    summary.to_csv(out/"variance.csv",index=False)
    selected=folds[folds.model=="boost_fusion"]
    selected.pivot(index="heat",columns="gas",values="r2")[list(GASES)].to_csv(out/"twenty_r2.csv")
    reference=pd.read_csv(RESULTS/"q3/curve_transfer/folds.csv").query("scenario=='H0_q2' and model=='tree_fusion'")
    comparison=selected.merge(reference,on=["heat","gas"],suffixes=("_boost","_tree"),validate="one_to_one")
    assert (comparison.n_boost==comparison.n_tree).all()
    assert (comparison.lead_boost==comparison.lead_tree).all()
    comparison["r2_gain"]=comparison.r2_boost-comparison.r2_tree
    comparison[["heat","gas","lead_boost","n_boost","r2_boost","r2_tree","r2_gain"]].to_csv(out/"comparison.csv",index=False)
    (out/"runtime.json").write_text(json.dumps(dict(python=platform.python_version(),sklearn=sklearn.__version__,
        seed=20261127,model="HistGradientBoostingRegressor",iterations=160,leaves=7,
        min_leaf=30,learning_rate=.04,l2=10,correction_weight=.5,early_stopping=False),indent=2),encoding="utf-8")
    experiment.clear()


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=RESULTS/"q3/boosted_curve/rerun")
    with threadpool_limits(1):
        run(parser.parse_args().output)
