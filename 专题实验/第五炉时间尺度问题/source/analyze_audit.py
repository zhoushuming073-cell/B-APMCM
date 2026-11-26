#!/usr/bin/env python3
"""Audit ten heats and compare 6C02057 under H0=0.5 s and H1=1.0 s."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
from scipy.ndimage import uniform_filter1d
from analysis_helpers import *

p=argparse.ArgumentParser(); p.add_argument('--cache-dir',default='cache'); p.add_argument('--out',default='results'); a=p.parse_args()
cache=Path(a.cache_dir).resolve(); out=Path(a.out).resolve(); tables=out/'tables'; features_dir=out/'features'
tables.mkdir(parents=True,exist_ok=True); features_dir.mkdir(parents=True,exist_ok=True)
D={}
for f in sorted(cache.glob('heats_*.npz')):
 z=np.load(f); h=f.stem.replace('heats_',''); D[h]={'seq':z['seq'],'X':z['X'],'proc':z['proc'],'wl':z['wl']}
if len(D)!=10: raise ValueError(f'Expected 10 caches, found {len(D)}')
TARGET='6C02057'; all_pos=np.concatenate([np.diff(d['proc'][:,0])[np.diff(d['proc'][:,0])>0] for d in D.values()]); quantum=float(np.median(all_pos))
rows=[]; updates=[]; deltarows=[]; qdiag=[]; flowrows=[]; gasdynrows=[]; gasintrows=[]; specrows=[]; assocrows=[]; lagrows=[]; crossrows=[]; cvrows=[]; spectral={}; target={}
for h,d in D.items():
 seq,X,P,wl=d['seq'],d['X'],d['proc'],d['wl']; O=P[:,0]; delta=np.diff(O); pos=delta[delta>0]; ix=np.flatnonzero(delta>0)+1; gaps=np.diff(np.r_[0,ix])
 rows.append({'heat':h,'rows':len(seq),'sequence_contiguous':bool(np.array_equal(seq,np.arange(1,len(seq)+1))),
  'duplicate_sequence':int(pd.Series(seq).duplicated().sum()),'missing_cells':int(np.isnan(np.c_[seq,X,P]).sum()),
  'duplicate_rows_excluding_sequence':int(pd.DataFrame(np.c_[X,P]).duplicated().sum()),
  'oxygen_start_m3':float(O[0]),'oxygen_end_m3':float(O[-1]),'oxygen_net_m3':float(O[-1]-O[0]),
  'duration_official_s':(len(O)-1)*.5,'mean_flow_official_m3_s':float((O[-1]-O[0])/((len(O)-1)*.5)),
  'flat_steps':int((delta==0).sum()),'positive_steps':int((delta>0).sum()),'negative_steps':int((delta<0).sum()),
  'flat_fraction':float(np.mean(delta==0)),'positive_increment_median_m3':float(np.median(pos)),
  'positive_increment_mean_m3':float(np.mean(pos)),'positive_increment_p90_m3':float(np.quantile(pos,.9)),
  'update_gap_median_rows':float(np.median(gaps)),'update_gap_p10_rows':float(np.quantile(gaps,.1)),
  'update_gap_p90_rows':float(np.quantile(gaps,.9)),'share_gap_1row':float(np.mean(gaps==1)),
  'share_gap_2rows':float(np.mean(gaps==2)),'share_gap_ge3rows':float(np.mean(gaps>=3)),
  'CO2_zero_rows':int((P[:,1]==0).sum()),'CO_zero_rows':int((P[:,2]==0).sum())})
 for k,ixrow in enumerate(ix):
  gap=int(gaps[k]); val=float(delta[ixrow-1]); updates.append({'heat':h,'row':int(ixrow+1),'delta_oxygen_m3':val,'gap_rows':gap,'interval_flow_H0_m3_s':val/(gap*.5)})
 for j,v in enumerate(delta): deltarows.append({'heat':h,'row':j+2,'delta_oxygen_m3':float(v)})
 # near-multiple-of-counter-increment proxy; descriptive only, not a confirmed meter quantum
 k=np.rint(pos/quantum).astype(int); err=np.abs(pos-k*quantum)
 qrows_heat={'heat':h,'n_positive_updates':int(len(pos)),'positive_delta_mean':float(np.mean(pos)),'positive_delta_median':float(np.median(pos)),'positive_delta_p90':float(np.quantile(pos,.9)),'reference_quantum_m3':quantum,'share_increment_2plus_quanta':float(np.mean(k>=2)),
  'share_increment_3plus_quanta':float(np.mean(k>=3)),'share_within_0p1_of_integer_quanta':float(np.mean(err<=.1)),
  'delta_gap_spearman':float(spearmanr(pos,gaps).statistic)}
 # gaps conditional on event increments
 for g in sorted(set(gaps)):
  x=pos[gaps==g]; qrows_heat[f'gap{g}_n']=len(x); qrows_heat[f'gap{g}_delta_median']=float(np.median(x))
 qdiag.append(qrows_heat)
 dt_opts={'H0':.5,'H1':1.0} if h==TARGET else {'official':.5}
 feat=spectral_features(X,wl); spectral[h]=feat
 pd.DataFrame({'sequence':seq,**feat}).to_csv(features_dir/f'{h}_spectral_features.csv',index=False)
 for clock,dt in dt_opts.items():
  q5=flow_slope(O,dt,5.,False); q10=flow_slope(O,dt,10.,False); qdirect=np.r_[delta/dt,np.nan]
  for name,q in [('direct_difference',qdirect),('centered_slope_5s',q5),('centered_slope_10s',q10)]:
   v=q[np.isfinite(q)]; flowrow={'heat':h,'clock':clock,'dt_s':dt,'method':name,'mean':float(v.mean()),
    'p10':float(np.quantile(v,.1)),'p50':float(np.quantile(v,.5)),'p90':float(np.quantile(v,.9)),
    'p95':float(np.quantile(v,.95)),'max':float(v.max()),'sd':float(v.std())}; flowrows.append(flowrow)
  for j,g in [(1,'CO2'),(2,'CO')]:
   y=P[:,j]; ac=acf(y,200); e=decay_lag(ac)
   gasdynrows.append({'heat':h,'clock':clock,'gas':g,'acf_e_fold_rows':e,'acf_e_fold_s':e*dt if np.isfinite(e) else np.nan,**psd_summary(y,dt),**event_duration(y,dt)})
   gasintrows.append({'heat':h,'clock':clock,'gas':g,'integral_m3':float(np.trapezoid(y,dx=dt))})
   ok=np.isfinite(q10)&np.isfinite(y); assocrows.append({'heat':h,'clock':clock,'gas':g,'flow_median_m3_s':float(np.nanmedian(q10)),
    'gas_median_m3_s':float(np.median(y)),'spearman_flow_gas':float(spearmanr(q10[ok],y[ok]).statistic),
    'pearson_flow_gas':float(pearsonr(q10[ok],y[ok]).statistic)})
   lagrows.append({'heat':h,'clock':clock,'gas':g,'flow_median_m3_s':float(np.nanmedian(q10)),**lag_estimate(q10,y,dt,20261126+sum(map(ord,h+clock+g)))})
  cO2=float(np.trapezoid(P[:,1]+.5*P[:,2],dx=dt)); total=float(np.trapezoid(P[:,1]+P[:,2],dx=dt))
  gasintrows.append({'heat':h,'clock':clock,'gas':'CO2_plus_CO','integral_m3':total,'carbon_o2_equivalent_m3':cO2,
   'oxygen_to_carbon_o2_equiv_ratio':float((O[-1]-O[0])/cO2) if cO2>0 else np.nan})
  for name in ['spectral_integral','spectral_p95','spectral_centroid','peak1_height','peak1_position_nm','peak1_width_nm','peak1_area']:
   if name not in feat: continue
   ac=acf(feat[name],200); e=decay_lag(ac)
   specrows.append({'heat':h,'clock':clock,'signal':name,'acf_e_fold_rows':e,'acf_e_fold_s':e*dt if np.isfinite(e) else np.nan,
    **psd_summary(feat[name],dt),**event_duration(feat[name],dt)})
  specrows.append({'heat':h,'clock':clock,'signal':'adjacent_cosine','adjacent_similarity_median':float(np.nanmedian(feat['adjacent_cosine']))})
  if h==TARGET: target[clock]={'q5':q5,'q10':q10,'qdirect':qdirect,'t':np.arange(len(O))*dt}
 if h==TARGET:
  for feature in ['spectral_integral','spectral_p95','spectral_centroid','peak1_height','peak1_position_nm','peak1_width_nm','peak1_area']:
   for gas,j in [('CO2',1),('CO',2)]:
    for clock,dt in [('H0',.5),('H1',1.)]:
     crossrows.append({'feature':feature,'gas':gas,'clock':clock,**lag_estimate(feat[feature],P[:,j],dt,20261126+len(crossrows))})
 print('audited',h,flush=True)

pd.DataFrame(rows).to_csv(tables/'oxygen_counter_audit.csv',index=False)
pd.DataFrame(updates).to_csv(tables/'oxygen_update_events.csv',index=False)
pd.DataFrame(deltarows).to_csv(tables/'oxygen_delta_rows.csv',index=False)
pd.DataFrame(qdiag).to_csv(tables/'counter_quantization_diagnostics.csv',index=False)
pd.DataFrame(flowrows).to_csv(tables/'oxygen_flow_summary.csv',index=False)
pd.DataFrame(gasdynrows).to_csv(tables/'gas_dynamics_summary.csv',index=False)
pd.DataFrame(gasintrows).to_csv(tables/'gas_integrals_and_balance.csv',index=False)
pd.DataFrame(assocrows).to_csv(tables/'flow_gas_associations.csv',index=False)
pd.DataFrame(lagrows).to_csv(tables/'lag_estimates.csv',index=False)
pd.DataFrame(specrows).to_csv(tables/'spectral_dynamics_summary.csv',index=False)
pd.DataFrame(crossrows).to_csv(tables/'spectral_gas_relative_lags.csv',index=False)
# detailed row-gap distribution
upd=[]
for h,d in D.items():
 delta=np.diff(d['proc'][:,0]); ix=np.flatnonzero(delta>0)+1; gaps=np.diff(np.r_[0,ix])
 for g in gaps: upd.append({'heat':h,'gap_rows':int(g),'gap_s_official':float(g*.5),'gap_s_H1_if_fifth':float(g) if h==TARGET else np.nan})
pd.DataFrame(upd).to_csv(tables/'oxygen_update_gap_events.csv',index=False)
# target nonnegative integrated B-spline CV under both clocks
cv=[]
for clock,dt in [('H0',.5),('H1',1.)]:
 _,_,r=fit_flow_cv(D[TARGET]['proc'][:,0],dt); r.update({'heat':TARGET,'clock':clock}); cv.append(r)
pd.DataFrame(cv).to_csv(tables/'fifth_oxygen_spline_cv.csv',index=False)
# pure row-based lag rescaling; identical filters and row lag, seconds are relabelled
z=D[TARGET]; O=z['proc'][:,0]; G=z['proc'][:,1:3]; qrow=flow_slope(O,1.,20.,False); pure=[]
for gas,j in [('CO2',0),('CO',1)]:
 qs=uniform_filter1d(np.nan_to_num(qrow,nan=np.nanmedian(qrow)),size=10,mode='nearest'); gs=uniform_filter1d(G[:,j],size=10,mode='nearest')
 cc=corr_lags(np.diff(qs,prepend=qs[0]),np.diff(gs,prepend=gs[0]),80); k=int(np.nanargmax(np.abs(cc)))
 pure.append({'gas':gas,'best_lag_rows':k,'corr':float(cc[k]),'H0_lag_s':k*.5,'H1_lag_s':k,'interpretation':'same row-index lag; physical seconds double by relabelling only'})
pd.DataFrame(pure).to_csv(tables/'pure_time_rescaling_lag_invariance.csv',index=False)
# summary
c=pd.DataFrame(rows).set_index('heat'); f=c.loc[TARGET]; o=c.drop(index=TARGET)
summary={'target':TARGET,'n_heats':10,'total_rows':int(c.rows.sum()),'all_sequences_contiguous':bool(c.sequence_contiguous.all()),
 'target_rows':int(f.rows),'target_start_m3':float(f.oxygen_start_m3),'target_end_m3':float(f.oxygen_end_m3),'target_net_m3':float(f.oxygen_net_m3),
 'target_duration_H0_s':(int(f.rows)-1)*.5,'target_duration_H1_s':(int(f.rows)-1),
 'target_mean_flow_H0_m3_s':float(f.oxygen_net_m3/((int(f.rows)-1)*.5)),
 'target_mean_flow_H1_m3_s':float(f.oxygen_net_m3/(int(f.rows)-1)),
 'target_positive_increment_median':float(f.positive_increment_median_m3),'target_positive_increment_mean':float(f.positive_increment_mean_m3),
 'target_update_gap_median_rows':float(f.update_gap_median_rows),
 'other_official_mean_flow_range':[float((o.oxygen_net_m3/((o.rows-1)*.5)).min()),float((o.oxygen_net_m3/((o.rows-1)*.5)).max())],
 'absolute_sample_interval_identifiable_from_rows_alone':False,
 'conclusion':'Retain H0 as official primary specification and H1 as sensitivity scenario; current attachment does not identify the fifth heat absolute clock.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
