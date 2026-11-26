"""Verify the imported experiment batch and recalculate its gas-assisted 5 s control."""
import hashlib
import io
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from python_calamine import CalamineWorkbook
from scipy.stats import spearmanr
from scipy.interpolate import BSpline
from scipy.optimize import lsq_linear
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT.parent/"APMCM_Full_Experiments"/"APMCM_Full_Experiments"
OUT=ROOT/"outputs/tables"
SHA="cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727"

def verified_series():
    archive=next(p for p in ROOT.parent.rglob("Data.zip") if file_hash(p)==SHA)
    checks=[];relations=[];data={}
    with ZipFile(archive) as z:
        for name in sorted(n for n in z.namelist() if n.lower().endswith(".xlsb")):
            h=Path(name).stem.split()[-1]
            blob=z.read(name)
            rows=CalamineWorkbook.from_filelike(io.BytesIO(blob)).get_sheet_by_index(0).iter_rows()
            header=next(rows);a=np.asarray(list(rows),float)
            seq=a[:,0];O=a[:,-3];G=a[:,-2:];t=(seq-1)*.5
            curve=pd.read_csv(PACKAGE/f"curves_{h}.csv")
            with np.load(PACKAGE/"cleaned_scenarios"/f"{h}_scenario.npz") as n:
                assert np.array_equal(a[:,1:-3],n["raw"]) and np.array_equal(seq,n["seq"])
                assert np.array_equal(np.array(header[1:-3],float),n["wl"])
                sat=n["raw"]==65535;near=np.zeros_like(sat)
                near[:,1:]|=sat[:,:-1];near[:,:-1]|=sat[:,1:]
                screen=n["screen"];repair=n["approved"];sr=n["saturation_repaired"];scenario=n["scenario"]
                assert np.array_equal(sat,n["saturated"]) and not (repair&~screen).any() and not (sr&~sat).any()
                assert np.array_equal(scenario[~(repair|sr)],n["raw"][~(repair|sr)])
                checks.append(dict(heat=h,sha256=hashlib.sha256(blob).hexdigest(),n=len(t),screen=screen.sum(),
                    approved=repair.sum(),sat_points=sat.sum(),sat_estimated=sr.sum(),sat_estimate_below_cap=(scenario[sr]<65535).sum(),
                    sat_horizontal_isolated=(sat&~near).sum(),gas_CO2_max=G[:,0].max(),gas_CO_max=G[:,1].max()))
            assert np.array_equal(t,curve.time_s) and np.allclose(curve[["oxygen_cumulative","CO2_flow","CO_flow"]],a[:,-3:],rtol=0,atol=1e-9)
            E=np.vstack([np.zeros(2),np.cumsum(.25*(G[1:]+G[:-1]),axis=0)])
            assert np.allclose(curve[["CO2_cumulative","CO_cumulative"]],E)
            for j,gas in enumerate(["CO2","CO"]):
                relations.append(dict(heat=h,gas=gas,volume=E[-1,j],rho_O_G=spearmanr(O,G[:,j]).statistic,
                    rho_O_E=spearmanr(O,E[:,j]).statistic,rho_dO_dG=spearmanr(np.diff(O),np.diff(G[:,j])).statistic,
                    rho_dO_dE=spearmanr(np.diff(O),np.diff(E[:,j])).statistic))
            data[h]=t,O,G
    assert len(data)==10
    pd.DataFrame(checks).to_csv(OUT/"full_batch_source_checks.csv",index=False)
    pd.DataFrame(relations).to_csv(OUT/"full_batch_oxygen_gas.csv",index=False)
    return data

def file_hash(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()

def refit_curves(data):
    params=pd.read_csv(PACKAGE/"q2_flow_nested_cv.csv").query("model == 'nonnegative_integral_spline'").set_index("heat")
    curves=[]
    for heat,(t,oxygen,gas) in data.items():
        spacing=params.loc[heat,"spacing_s"];penalty=params.loc[heat,"penalty"]
        knots=np.r_[[t[0]]*4,np.arange(t[0]+spacing,t[-1],spacing),[t[-1]]*4]
        count=len(knots)-4;spline=BSpline(knots,np.eye(count),3)
        basis=spline(t);anti=spline.antiderivative();integrated=anti(t)-anti(t[0])
        regularizer=np.sqrt(penalty)*np.diff(np.eye(count),n=2,axis=0)
        weight=np.ones(len(t))
        for _ in range(2):
            design=np.vstack([integrated*np.sqrt(weight[:,None]),regularizer])
            target=np.r_[(oxygen-oxygen[0])*np.sqrt(weight),np.zeros(len(regularizer))]
            result=lsq_linear(design,target,bounds=(0,np.inf),tol=1e-7,max_iter=100)
            if not result.success: raise RuntimeError(f"Flow refit failed: {heat}")
            fitted=oxygen[0]+integrated@result.x
            residual=oxygen-fitted
            scale=max(.1,1.4826*np.median(abs(residual-np.median(residual))))
            weight=np.minimum(1,1.345*scale/np.maximum(abs(residual),1e-9))
        flow=basis@result.x
        assert min(flow)>=-1e-8 and min(np.diff(fitted))>=-1e-7
        emission=np.vstack([np.zeros(2),np.cumsum(.25*(gas[1:]+gas[:-1]),axis=0)])
        curves.append(pd.DataFrame(dict(heat=heat,time_s=t,oxygen_measured_m3=oxygen,oxygen_fit_m3=fitted,
            oxygen_flow_m3s=flow,CO2_flow_m3s=gas[:,0],CO_flow_m3s=gas[:,1],
            CO2_cumulative_m3=emission[:,0],CO_cumulative_m3=emission[:,1])))
    pd.concat(curves,ignore_index=True).to_csv(OUT/"full_batch_curves_refit.csv",index=False)

def features(t,O,G):
    # Equal grid/issue times to the imported experiment. No whole-heat imputation.
    q=np.zeros(len(O))
    for i in range(2,len(O)):
        j=max(0,i-19);u=t[j:i+1]-t[j:i+1].mean()
        q[i]=max(0,float(u@O[j:i+1]/(u@u)))
    t,O,G,q=t[::2],O[::2],G[::2],q[::2]
    issue=np.arange(40,len(t)-5)
    progress=np.column_stack([t[issue],O[issue]])
    last=G[issue]
    base=np.column_stack([last,progress])
    lags=np.column_stack([q[np.maximum(0,issue-k)] for k in [0,5,10,20,30,40]])
    histories=np.column_stack([[q[max(0,i-older):max(0,i-younger)].mean() if i>younger else q[i]
                                for i in issue] for younger,older in [(0,5),(5,15),(15,30),(30,40),(40,60)]])
    return dict(last=last,gas_ridge=last,progress=base,current=np.column_stack([base,lags[:,0]]),
                history=np.column_stack([base,histories]),laggrid=np.column_stack([base,lags])),G[issue+5]

def run_controls(data):
    f={h:features(*x) for h,x in data.items()};results=[]
    for test,(feat,y) in f.items():
        train=[h for h in f if h!=test]
        for j,gas in enumerate(["CO2","CO"]):
            for key in feat:
                alpha=np.nan
                if key=="last":pred=feat[key][:,j]
                else:
                    def predict(ids,v,hp):
                        X=np.vstack([f[h][0][key] for h in ids]);target=np.concatenate([f[h][1][:,j] for h in ids])
                        reg=make_pipeline(StandardScaler(),Ridge(alpha=hp)).fit(X,target)
                        return np.maximum(0,reg.predict(f[v][0][key]))
                    alphas=[.1,1,10,100]
                    losses=[np.mean([np.mean(abs(predict([h for h in train if h!=v],v,hp)-f[v][1][:,j]))
                                     for v in train]) for hp in alphas]
                    alpha=alphas[np.argmin(losses)];pred=predict(train,test,alpha)
                e=pred-y[:,j];results.append(dict(heat=test,gas=gas,model=key,alpha=alpha,n=len(e),
                    mae=np.abs(e).mean(),rmse=np.sqrt((e**2).mean())))
    df=pd.DataFrame(results);df.to_csv(OUT/"full_batch_causal_5s.csv",index=False)
    print(df.groupby(["gas","model"])[["mae","rmse"]].mean().round(5))

if __name__=="__main__":
    with threadpool_limits(1):
        data=verified_series()
        refit_curves(data)
        run_controls(data)
