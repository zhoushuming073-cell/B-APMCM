from __future__ import annotations
from pathlib import Path
import json
import numpy as np,pandas as pd

def run(results_dir,archive_dir,flow_profile_dir=None,round3_curves_path=None):
    r=Path(results_dir);a=Path(archive_dir)
    # Compare the archived compact 24-row second-round reference against the regenerated 40-test family.
    old=pd.read_csv(a/'key_actual_results.csv');new=pd.read_csv(r/'round2_repaired_40_test_family.csv')
    keys=['scenario','method','gas','window_s','offset_s']
    m=old.merge(new,on=keys,suffixes=('_arch','_new'),how='left',indicator=True)
    fields=[('median_lag_s_arch','median_lag_s_new'),('high_minus_low_lag_s','median_high_minus_low_lag_s'),
            ('n_heats_contrast_arch','n_heats_contrast_new'),('n_heats_high_flow_shorter_arch','n_heats_high_flow_shorter_new'),
            ('sign_p','sign_test_one_sided_p'),('holm_p','holm_p_40')]
    rows=[]
    for _,x in m.iterrows():
        checks=[];not_archived=0
        for oldc,newc in fields:
            ov=pd.to_numeric(pd.Series([x.get(oldc,np.nan)]),errors='coerce').iloc[0]
            nv=pd.to_numeric(pd.Series([x.get(newc,np.nan)]),errors='coerce').iloc[0]
            if pd.isna(ov):
                not_archived+=1;continue
            checks.append(bool(pd.notna(nv) and np.isclose(ov,nv,atol=1e-10)))
        rows.append(dict(scenario=x.scenario,method=x.method,gas=x.gas,window_s=x.window_s,offset_s=x.offset_s,
            archive_record_found=x['_merge']=='both',available_fields_compared=len(checks),available_fields_match=sum(checks),
            archive_fields_not_recorded=not_archived,all_available_fields_match=all(checks)))
    comp=pd.DataFrame(rows);comp.to_csv(r/'round2_archive_comparison.csv',index=False)
    # Compare the formal lag summary against the official frozen summary.
    oldlag=pd.read_csv(a/'lag_summary.csv');newlag=pd.read_csv(r/'actual_lag_summary.csv').query("metric=='pearson'")
    lm=oldlag.merge(newlag,on=['metric','scenario','window_s','gas'],suffixes=('_arch','_new'),how='left')
    lagfields=[('n_windows_arch','n_windows_new'),('n_heats_arch','n_heats_new'),('median_lag_s','lag_median_s'),
               ('lag_q25_s_arch','lag_q25_s_new'),('lag_q75_s_arch','lag_q75_s_new'),('supported_windows_q_le_0_1','bh_q_le_0_1_n')]
    lrows=[]
    for _,x in lm.iterrows():
        vals=[]
        for oc,nc in lagfields:
            ov=pd.to_numeric(pd.Series([x[oc]]),errors='coerce').iloc[0];nv=pd.to_numeric(pd.Series([x[nc]]),errors='coerce').iloc[0]
            vals.append(bool((pd.isna(ov) and pd.isna(nv)) or (pd.notna(ov) and pd.notna(nv) and np.isclose(ov,nv,atol=1e-10))))
        lrows.append(dict(scenario=x.scenario,window_s=x.window_s,gas=x.gas,field_match_count=sum(vals),fields_compared=len(vals),all_fields_match=all(vals)))
    pd.DataFrame(lrows).to_csv(r/'lag_summary_archive_comparison.csv',index=False)
    # Compare flow-bin lag medians (the official table does not contain the actual bin boundaries).
    oldflow=pd.read_csv(a/'flow_delay_summary.csv');bins=pd.read_csv(r/'actual_flow_bins.csv').query("metric=='pearson'")
    frows=[]
    for _,x in oldflow.iterrows():
        z=bins[(bins.scenario==x.scenario)&(bins.window_s==90)&(bins.gas==x.gas)].set_index('flow_bin')
        vals=[float(z.loc[k,'lag_median_s']) for k in ['low','middle','high']]
        refs=[float(x.low_flow_lag_median_s),float(x.middle_flow_lag_median_s),float(x.high_flow_lag_median_s)]
        frows.append(dict(scenario=x.scenario,gas=x.gas,rebuilt_low=vals[0],archived_low=refs[0],rebuilt_middle=vals[1],archived_middle=refs[1],rebuilt_high=vals[2],archived_high=refs[2],all_match=bool(np.allclose(vals,refs))))
    pd.DataFrame(frows).to_csv(r/'flow_delay_archive_comparison.csv',index=False)
    # Historical 180-s prewhitening version discrepancy.
    pd.read_csv(a/'prewhitened_180_version.csv').to_csv(r/'prewhitened_version_comparison.csv',index=False)
    if round3_curves_path:
        current_path=r.parent/'curves'/'final_curves.csv';ref_path=Path(round3_curves_path)
        cur=pd.read_csv(current_path);ref=pd.read_csv(ref_path)
        keys=['scenario','heat','row'];cm=ref.merge(cur,on=keys,suffixes=('_round3','_current'),how='outer',indicator=True,validate='one_to_one')
        numeric=[c for c in ref.columns if c not in keys and pd.api.types.is_numeric_dtype(ref[c])]
        dr=[]
        for c in numeric:
            delta=(pd.to_numeric(cm[c+'_current'],errors='coerce')-pd.to_numeric(cm[c+'_round3'],errors='coerce')).abs()
            dr.append(dict(column=c,max_abs_difference=float(delta.max()),mean_abs_difference=float(delta.mean()),n_nonzero=int((delta>0).sum())))
        pd.DataFrame(dr).to_csv(r/'round3_curve_numeric_difference.csv',index=False)
        current_sha=__import__('hashlib').sha256(current_path.read_bytes()).hexdigest();round3_sha=__import__('hashlib').sha256(ref_path.read_bytes()).hexdigest()
        pd.DataFrame([dict(current_curve_sha256=current_sha,round3_curve_sha256=round3_sha,rows_current=len(cur),rows_round3=len(ref),
            unique_key_current=not cur.duplicated(keys).any(),unique_key_round3=not ref.duplicated(keys).any(),
            unmatched_rows=int((cm['_merge']!='both').sum()),max_abs_numeric_difference=float(max(x['max_abs_difference'] for x in dr)),
            mean_abs_numeric_difference=float(np.mean([x['mean_abs_difference'] for x in dr])),
            current_matches_relationship_study_hash=current_sha=='9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49',
            round3_reference_hash_matches_round3_manifest=round3_sha=='c89a8e34dc5f83385d6c06d356d2821f53e67b63197aa00b6dca9cf0f837f306')]).to_csv(r/'curve_reconstruction_archive_comparison.csv',index=False)
        qa=json.loads((r/'qa_manifest.json').read_text()) if (r/'qa_manifest.json').exists() else {}
        pd.DataFrame([
          dict(item='Data.zip SHA256',expected='cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727',observed=qa.get('source_data_zip_sha256','NOT_AVAILABLE'),status='MATCH' if qa.get('source_data_zip_sha256')=='cb662ce5ea17c247d0733ec8a4245fec37abf82d823a0ddc08597b300d95a727' else 'MISMATCH'),
          dict(item='process.csv GitHub blob SHA1',expected='25a8f164274ff96e868ff447d28a4d1c80353b50',observed=qa.get('process_csv_github_blob_sha1','NOT_AVAILABLE'),status='MATCH' if qa.get('process_csv_github_blob_sha1')=='25a8f164274ff96e868ff447d28a4d1c80353b50' else 'MISMATCH'),
          dict(item='Q2.1 manifest legacy process CSV SHA256',expected='f82e726b03f8eceb148a25440f34fe8385af28e95a7f551149b513ecf560c789',observed=qa.get('process_sha256','NOT_AVAILABLE'),status='MISMATCH_RETAINED_FOR_AUDIT'),
          dict(item='Relationship study frozen curve SHA256',expected='9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49',observed=current_sha,status='MATCH' if current_sha=='9bdc1fe5d962e5ed3b2c7ae513aade8094ea2b3e95a7f405a12c485d6099fd49' else 'MISMATCH'),
          dict(item='Round3 input curve SHA256',expected='c89a8e34dc5f83385d6c06d356d2821f53e67b63197aa00b6dca9cf0f837f306',observed=round3_sha,status='MATCH' if round3_sha=='c89a8e34dc5f83385d6c06d356d2821f53e67b63197aa00b6dca9cf0f837f306' else 'MISMATCH'),
          dict(item='Round4 input curve SHA256',expected='f7280a7ee7276c6281b2c5f85b9fb003e73ea51ad8473656873a95f07fd016b4',observed='NOT_AVAILABLE',status='HISTORICAL_INPUT_NOT_AVAILABLE')]).to_csv(r/'historical_input_hash_audit.csv',index=False)
    if flow_profile_dir:
        fp=Path(flow_profile_dir);official=pd.read_csv(a/'flow_summary.csv');rebuilt=pd.read_csv(fp/'flow_summary.csv')
        fm=official.merge(rebuilt,on=['scenario','heat'],suffixes=('_official','_frozen'),how='left')
        numeric=[c for c in official.columns if c not in ['scenario','heat','fit_success']]
        frows=[]
        for _,x in fm.iterrows():
            checks=[bool(pd.notna(x.get(c+'_frozen')) and np.isclose(float(x[c+'_official']),float(x[c+'_frozen']),rtol=1e-10,atol=1e-10)) for c in numeric]
            frows.append(dict(scenario=x.scenario,heat=x.heat,numeric_fields_match=sum(checks),numeric_fields_compared=len(checks),
                fit_success_match=bool(x.fit_success_official==x.fit_success_frozen),all_fields_match=all(checks) and bool(x.fit_success_official==x.fit_success_frozen)))
        pd.DataFrame(frows).to_csv(r/'flow_profile_archive_comparison.csv',index=False)
    return comp,lm

if __name__=='__main__':
    import sys;run(sys.argv[1],sys.argv[2],sys.argv[3] if len(sys.argv)>3 else None,sys.argv[4] if len(sys.argv)>4 else None)
