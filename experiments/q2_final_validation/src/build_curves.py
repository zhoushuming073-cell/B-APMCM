from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd
from scipy.interpolate import PchipInterpolator
from flow_v4 import run as run_flow

def fit_k2(t,y):
    t=np.asarray(t,float);y=np.asarray(y,float)
    knots=np.unique(np.r_[t[0],np.arange(t[0]+2.,t[-1],2.),t[-1]])
    centers=[t[0]];means=[y[0]]
    for knot in knots[1:-1]:
        take=np.abs(t-knot)<=1.
        if take.sum()>=2: centers.append(float(t[take].mean()));means.append(float(y[take].mean()))
    centers.append(t[-1]);means.append(y[-1])
    return PchipInterpolator(centers,np.maximum(means,0),extrapolate=False)

def build(input_csv,outdir,archive_path=None,flow_profile_path=None):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True); flowdir=outdir/'flow'
    if flow_profile_path:
        fp_path=Path(flow_profile_path);profile=pd.read_csv(fp_path);flowdir=fp_path.parent
        qa_path=flowdir/'flow_validation_checks.csv'
        if not qa_path.exists():raise FileNotFoundError(f'Missing frozen Q2.1 validation checks: {qa_path}')
        qa=pd.read_csv(qa_path)
    else:
        profile,_,qa=run_flow(Path(input_csv),flowdir)
    src=pd.read_csv(input_csv).sort_values(['attachment','row']); curves=[]
    for scenario in ('H0','H1'):
      for heat,d in src.groupby('heat',sort=False):
        d=d.sort_values('row'); dt=1. if scenario=='H1' and heat=='6C02057' else .5; t=(d.row.to_numpy(float)-1)*dt
        fp=profile[(profile.scenario==scenario)&(profile.heat==heat)].sort_values('row')
        vals={}
        for gas in ('CO2','CO'):
            p=fit_k2(t,d[gas].to_numpy(float)); yhat=p(t); anti=p.antiderivative();cum=anti(t)-anti(t[0]);vals[gas]=(yhat,cum)
        curves.append(pd.DataFrame({'scenario':scenario,'heat':heat,'row':d.row.to_numpy(int),'time_s':t,
          'oxygen_cumulative_raw_m3':d.O.to_numpy(float),'oxygen_cumulative_fit_m3':fp.cumulative_O_fit_m3.to_numpy(float),
          'oxygen_flow_analytic_m3_s':fp.flow_analytic_best_m3_s.to_numpy(float),
          'CO2_flow_raw_m3_s':d.CO2.to_numpy(float),'CO2_flow_k2_m3_s':vals['CO2'][0],'CO2_cumulative_k2_m3':vals['CO2'][1],
          'CO_flow_raw_m3_s':d.CO.to_numpy(float),'CO_flow_k2_m3_s':vals['CO'][0],'CO_cumulative_k2_m3':vals['CO'][1]}))
    curves=pd.concat(curves,ignore_index=True); path=outdir/'final_curves.csv';curves.to_csv(path,index=False)
    m={'source_process_csv_sha256':hashlib.sha256(Path(input_csv).read_bytes()).hexdigest(),'source_archive_sha256':hashlib.sha256(Path(archive_path if archive_path else '/mnt/data/Data.zip').read_bytes()).hexdigest(),'rows':len(curves),'heats':int(curves.heat.nunique()),'scenarios':sorted(curves.scenario.unique()),'curve_csv_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'flow_validation_pass':bool(qa[['pass_nonnegative','pass_monotonic','pass_mass_balance_0_1m3']].all().all()),'flow_validation_rows':len(qa),'flow_profile_sha256':hashlib.sha256((Path(flow_profile_path) if flow_profile_path else flowdir/'flow_profile.csv').read_bytes()).hexdigest(),'flow_profile_source':'frozen_q2.1_profile' if flow_profile_path else 'recomputed_from_process_csv'}
    (outdir/'manifest.json').write_text(json.dumps(m,indent=2))
    return curves,m,qa
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input',default='/mnt/data/q23_final_validation/inputs/process.csv');ap.add_argument('--out',default='/mnt/data/q23_final_validation/results/rebuild');ap.add_argument('--flow-profile',default=None);a=ap.parse_args()
 c,m,q=build(a.input,a.out,flow_profile_path=a.flow_profile);print(m);print(q.groupby('scenario')[['pass_nonnegative','pass_monotonic','pass_mass_balance_0_1m3']].all())
