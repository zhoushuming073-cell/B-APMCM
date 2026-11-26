from pathlib import Path
import json, math, zipfile, shutil, hashlib, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, correlate
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
from sklearn.isotonic import IsotonicRegression

ROOT = Path(__file__).resolve().parents[1]
CURVE_DIR = ROOT/'inputs'
SPEC_FILE = ROOT/'inputs/harmonized_spectra.npz'
OUT = ROOT/'rerun_output'
if OUT.exists():
    shutil.rmtree(OUT)
for name in ['tables','figures','code','inputs']:
    (OUT/name).mkdir(parents=True, exist_ok=True)

HEATS = ['6B03205','6B03206','6B03270','6B03271','6C02057','6C02058','6C02059','6C02060','6C02107','6C02108']
GASES = ['CO2','CO']
BASE_DT = 0.5
MAX_LAG_S = 60.0
WINDOW_S = 120.0
STEP_S = 60.0
BOOT_REPS = 60
NULL_REPS = 199
SEED = 20261126

curves = {}
for h in HEATS:
    p = CURVE_DIR/f'curves_{h}.csv'
    if not p.exists():
        raise FileNotFoundError(p)
    curves[h] = pd.read_csv(p)

# Helpers
def rolling_slope(y, dt, window_s=10.0, centered=True):
    y = np.asarray(y, float)
    n = len(y)
    k = max(5, int(round(window_s/dt)))
    if k % 2 == 0:
        k += 1
    out = np.full(n, np.nan)
    half = k//2
    for i in range(n):
        if centered:
            a, b = max(0, i-half), min(n, i+half+1)
        else:
            a, b = max(0, i-k+1), i+1
        if b-a < 4:
            continue
        xx = np.arange(a,b,dtype=float)*dt
        yy = y[a:b]
        xc = xx-xx.mean()
        den = np.dot(xc,xc)
        if den > 0:
            out[i] = max(0.0, float(np.dot(xc,yy-yy.mean())/den))
    # Fill edge gaps only by nearest valid estimate; do not extrapolate trend.
    good = np.isfinite(out)
    if good.any():
        out[:np.where(good)[0][0]] = out[good][0]
        out[np.where(good)[0][-1]+1:] = out[good][-1]
        out = pd.Series(out).interpolate(limit_direction='both').to_numpy()
    return out

def odd_window(seconds, dt, min_n=5):
    w = max(min_n, int(round(seconds/dt)))
    if w % 2 == 0:
        w += 1
    return w

def smooth_derivative(y, dt, smooth_s=5.0):
    y = np.asarray(y, float)
    w = odd_window(smooth_s, dt)
    w = min(w, len(y) if len(y)%2 else len(y)-1)
    if w < 5:
        return np.gradient(y, dt)
    ys = savgol_filter(y, window_length=w, polyorder=2, mode='interp')
    return np.gradient(ys, dt)

def corr_lags(x, y, dt, maxlag_s=MAX_LAG_S, min_overlap_frac=0.5):
    x = np.asarray(x,float); y=np.asarray(y,float)
    n=min(len(x),len(y)); x=x[:n]; y=y[:n]
    maxk=min(int(round(maxlag_s/dt)), n-5)
    min_n=max(40, int(math.ceil(n*min_overlap_frac)))
    # FFT cross-correlation plus prefix sums gives all lag-wise Pearson r values
    # in O(n log n), rather than repeatedly scanning the full series.
    cross=correlate(y,x,mode='full',method='fft')[n-1:n+maxk]
    px=np.r_[0.0,np.cumsum(x)]; py=np.r_[0.0,np.cumsum(y)]
    px2=np.r_[0.0,np.cumsum(x*x)]; py2=np.r_[0.0,np.cumsum(y*y)]
    ks=np.arange(maxk+1); m=n-ks
    sx=px[n-ks]; sy=py[n]-py[ks]
    sx2=px2[n-ks]; sy2=py2[n]-py2[ks]
    cov=cross-sx*sy/m
    vx=np.maximum(0.0,sx2-sx*sx/m); vy=np.maximum(0.0,sy2-sy*sy/m)
    den=np.sqrt(vx*vy)
    rs=np.divide(cov,den,out=np.full_like(cov,np.nan,dtype=float),where=(den>1e-12)&(m>=min_n))
    rs=np.clip(rs,-1,1)
    return ks*dt,rs

def peak_summary(lags, rs):
    ok=np.isfinite(rs)
    if not ok.any():
        return dict(best_abs_lag_s=np.nan,best_abs_r=np.nan,best_abs_sign=0,
                    best_positive_lag_s=np.nan,best_positive_r=np.nan,
                    best_negative_lag_s=np.nan,best_negative_r=np.nan,
                    zero_lag_r=np.nan,peak_width90_s=np.nan,boundary=True)
    ll=lags[ok]; rr=rs[ok]
    ia=int(np.argmax(np.abs(rr))); ip=int(np.argmax(rr)); ine=int(np.argmin(rr))
    peak=abs(rr[ia]); thr=.9*peak
    near=np.where(np.abs(rr)>=thr)[0]
    # width of contiguous neighborhood around max abs peak
    left=ia; right=ia
    while left-1>=0 and abs(rr[left-1])>=thr: left-=1
    while right+1<len(rr) and abs(rr[right+1])>=thr: right+=1
    return dict(best_abs_lag_s=float(ll[ia]),best_abs_r=float(rr[ia]),best_abs_sign=int(np.sign(rr[ia])),
                best_positive_lag_s=float(ll[ip]),best_positive_r=float(rr[ip]),
                best_negative_lag_s=float(ll[ine]),best_negative_r=float(rr[ine]),
                zero_lag_r=float(rr[0]),peak_width90_s=float(ll[right]-ll[left]),
                boundary=bool(ia==0 or ia==len(rr)-1))

def circular_shift_null(x,y,dt,maxlag_s=MAX_LAG_S,reps=NULL_REPS,seed=0):
    rng=np.random.default_rng(seed)
    n=min(len(x),len(y)); x=np.asarray(x[:n],float); y=np.asarray(y[:n],float)
    lags,rs=corr_lags(x,y,dt,maxlag_s)
    obs=np.nanmax(np.abs(rs)) if np.isfinite(rs).any() else np.nan
    min_shift=min(n//3, int(round((maxlag_s+30)/dt)))
    max_shift=n-min_shift
    if not np.isfinite(obs) or max_shift<=min_shift:
        return np.nan, np.nan, np.array([])
    null=[]
    for _ in range(reps):
        s=int(rng.integers(min_shift,max_shift))
        yr=np.roll(y,s)
        _,rn=corr_lags(x,yr,dt,maxlag_s)
        null.append(np.nanmax(np.abs(rn)) if np.isfinite(rn).any() else np.nan)
    null=np.asarray(null)
    p=(1+np.sum(null>=obs))/(1+np.sum(np.isfinite(null)))
    return float(obs),float(p),null

def paired_block_bootstrap(x,y,dt,block_s=120.0,reps=BOOT_REPS,seed=0,maxlag_s=MAX_LAG_S):
    rng=np.random.default_rng(seed)
    x=np.asarray(x,float); y=np.asarray(y,float)
    n=min(len(x),len(y)); x=x[:n]; y=y[:n]
    bl=max(int(round(block_s/dt)), int(round((maxlag_s+10)/dt)))
    bl=min(bl,n)
    nblocks=int(math.ceil(n/bl))
    lags_out=[]; r_out=[]; signs=[]
    for _ in range(reps):
        starts=rng.integers(0,n,size=nblocks)
        xb=[]; yb=[]
        for st in starts:
            ix=(st+np.arange(bl))%n
            xb.append(x[ix]); yb.append(y[ix])
        xx=np.concatenate(xb)[:n]; yy=np.concatenate(yb)[:n]
        ll,rr=corr_lags(xx,yy,dt,maxlag_s)
        if not np.isfinite(rr).any(): continue
        k=int(np.nanargmax(np.abs(rr)))
        lags_out.append(ll[k]); r_out.append(rr[k]); signs.append(np.sign(rr[k]))
    if not lags_out:
        return np.nan,np.nan,np.nan,np.nan,np.nan
    lags_out=np.asarray(lags_out); r_out=np.asarray(r_out); signs=np.asarray(signs)
    return (float(np.quantile(lags_out,.025)),float(np.quantile(lags_out,.975)),
            float(np.median(lags_out)),float(np.quantile(np.abs(r_out),.5)),
            float(np.mean(signs==np.sign(np.median(signs)))))

def bh_adjust(df, pcol='circular_shift_p'):
    df=df.copy(); valid=df[pcol].notna()
    df['BH_q']=np.nan
    if valid.sum():
        df.loc[valid,'BH_q']=multipletests(df.loc[valid,pcol].values,method='fdr_bh')[1]
    return df

# Construct scenario-specific signals and full-series lag diagnostics.
full_rows=[]; window_rows=[]; series_cache={}
for scenario in ['H0','H1']:
    for hi,h in enumerate(HEATS):
        print('processing',scenario,h,flush=True)
        d=curves[h]
        dt=1.0 if (scenario=='H1' and h=='6C02057') else BASE_DT
        O=d['oxygen_cumulative'].to_numpy(float)
        q=rolling_slope(O,dt,10.0,centered=True)
        qt=rolling_slope(O,dt,10.0,centered=False)
        for gas in GASES:
            g=d[f'{gas}_flow'].to_numpy(float)
            # Centered offline reference and trailing causal-reference sensitivity.
            x=smooth_derivative(q,dt,5.0)
            y=smooth_derivative(g,dt,5.0)
            xt=smooth_derivative(qt,dt,5.0)
            # Remove edge region so smoothing and lag overlap do not dominate.
            trim=min(int(round(20/dt)), max(0,(len(x)-100)//4))
            if trim:
                xx=x[trim:-trim]; yy=y[trim:-trim]; xct=xt[trim:-trim]
            else:
                xx=x; yy=y; xct=xt
            ll,rr=corr_lags(xx,yy,dt,MAX_LAG_S)
            sm=peak_summary(ll,rr)
            _,p,null=circular_shift_null(xx,yy,dt,MAX_LAG_S,NULL_REPS,SEED+hi*13+(0 if gas=='CO2' else 7)+(0 if scenario=='H0' else 1000))
            ci_lo,ci_hi,boot_med,boot_absr,boot_sign=paired_block_bootstrap(xx,yy,dt,120.0,BOOT_REPS,SEED+hi*17+(1 if gas=='CO' else 0)+(1000 if scenario=='H1' else 0))
            llt,rrt=corr_lags(xct,yy,dt,MAX_LAG_S)
            smt=peak_summary(llt,rrt)
            row={'scenario':scenario,'heat':h,'gas':gas,'dt_s':dt,'n_rows':len(d),
                 'q_median_m3s':float(np.nanmedian(q)),'q_p10_m3s':float(np.nanquantile(q,.1)),
                 'q_p90_m3s':float(np.nanquantile(q,.9)), 'gas_sd_m3s':float(np.std(g)),
                 **sm,'bootstrap_lag_ci_low_s':ci_lo,'bootstrap_lag_ci_high_s':ci_hi,
                 'bootstrap_lag_median_s':boot_med,'bootstrap_abs_r_median':boot_absr,
                 'bootstrap_sign_stability':boot_sign,
                 'circular_shift_p':p,'causal_q_best_abs_lag_s':smt['best_abs_lag_s'],
                 'causal_q_best_abs_r':smt['best_abs_r'],
                 'causal_q_boundary':smt['boundary']}
            full_rows.append(row)
            series_cache[(scenario,h,gas)]={'dt':dt,'q':q,'qt':qt,'x':x,'y':y,'g':g,'lags':ll,'corr':rr}
            # 120-second windows, 60-second stride, maximum local lag 40 seconds.
            wlen=max(80,int(round(WINDOW_S/dt))); step=max(40,int(round(STEP_S/dt)))
            maxlag_local=40.0
            for start in range(0,len(d)-wlen+1,step):
                end=start+wlen
                xw=x[start:end]; yw=y[start:end]; qw=q[start:end]
                l2,r2=corr_lags(xw,yw,dt,maxlag_local,min_overlap_frac=.55)
                ps=peak_summary(l2,r2)
                qsd=float(np.std(qw)); qmed=float(np.median(qw))
                informative=(qsd>=.15 and np.isfinite(ps['best_abs_r']) and abs(ps['best_abs_r'])>=.25 and not ps['boundary'])
                window_rows.append({'scenario':scenario,'heat':h,'gas':gas,
                    'window_start_s':float(start*dt),'window_end_s':float((end-1)*dt),
                    'flow_median_m3s':qmed,'flow_sd_m3s':qsd,
                    'best_lag_s':ps['best_abs_lag_s'],'best_corr':ps['best_abs_r'],
                    'best_sign':ps['best_abs_sign'],'peak_width90_s':ps['peak_width90_s'],
                    'boundary':ps['boundary'],'informative':bool(informative)})

full=pd.DataFrame(full_rows)
# Adjust within each scenario across all heat-gas tests.
full=pd.concat([bh_adjust(g) for _,g in full.groupby('scenario')],ignore_index=True)
full.to_csv(OUT/'tables/full_series_lag_diagnostics.csv',index=False)

windows=pd.DataFrame(window_rows)
# Define common flow support from nine non-fifth heats under official H0.
ref=windows[(windows.scenario=='H0')&(windows.heat!='6C02057')]
flow_lo=float(ref.flow_median_m3s.quantile(.01)); flow_hi=float(ref.flow_median_m3s.quantile(.99))
qcuts=ref.flow_median_m3s.quantile([1/3,2/3]).to_numpy()
windows['in_reference_flow_support']=windows.flow_median_m3s.between(flow_lo,flow_hi)
windows['flow_bin']=pd.cut(windows.flow_median_m3s,[-np.inf,qcuts[0],qcuts[1],np.inf],labels=['low','middle','high'])
windows.to_csv(OUT/'tables/local_window_lag_diagnostics.csv',index=False)

# Summaries by flow bin; report both all windows and informative windows.
bin_rows=[]
for scenario in ['H0','H1']:
    for gas in GASES:
        sub=windows[(windows.scenario==scenario)&(windows.gas==gas)]
        for subset_name,ss in [('all',sub),('informative',sub[sub.informative])]:
            for b in ['low','middle','high']:
                z=ss[ss.flow_bin==b]
                bin_rows.append({'scenario':scenario,'gas':gas,'window_set':subset_name,'flow_bin':b,
                    'n_windows':len(z),'n_heats':z.heat.nunique(),
                    'flow_median_m3s':z.flow_median_m3s.median() if len(z) else np.nan,
                    'lag_median_s':z.best_lag_s.median() if len(z) else np.nan,
                    'lag_q25_s':z.best_lag_s.quantile(.25) if len(z) else np.nan,
                    'lag_q75_s':z.best_lag_s.quantile(.75) if len(z) else np.nan,
                    'median_abs_corr':z.best_corr.median() if len(z) else np.nan,
                    'negative_corr_share':(z.best_sign<0).mean() if len(z) else np.nan})
bins=pd.DataFrame(bin_rows)
bins.to_csv(OUT/'tables/flow_bin_lag_summary.csv',index=False)

# Within-heat flow-lag slope: lag = heat intercept + beta * (flow - heat mean).
def within_slope(df):
    d=df.dropna(subset=['flow_median_m3s','best_lag_s']).copy()
    if len(d)<5 or d.heat.nunique()<3: return np.nan
    dx=d.flow_median_m3s-d.groupby('heat').flow_median_m3s.transform('mean')
    dy=d.best_lag_s-d.groupby('heat').best_lag_s.transform('mean')
    den=float(np.dot(dx,dx))
    return float(np.dot(dx,dy)/den) if den>0 else np.nan

slope_rows=[]; slope_boot_rows=[]
for scenario in ['H0','H1']:
    for gas in GASES:
        z=windows[(windows.scenario==scenario)&(windows.gas==gas)&(windows.informative)&(windows.in_reference_flow_support)].copy()
        # H0 fifth is outside common support in most cases; explicitly exclude it from generalizable slope.
        if scenario=='H0': z=z[z.heat!='6C02057']
        per=[]
        for h,g in z.groupby('heat'):
            if len(g)>=3:
                rho=spearmanr(g.flow_median_m3s,g.best_lag_s).statistic if g.flow_median_m3s.nunique()>1 and g.best_lag_s.nunique()>1 else np.nan
                per.append({'heat':h,'n_windows':len(g),'rho':rho,'slope_s_per_m3s':np.polyfit(g.flow_median_m3s,g.best_lag_s,1)[0] if g.flow_median_m3s.nunique()>1 else np.nan})
        per_df=pd.DataFrame(per)
        beta=within_slope(z)
        # Cluster bootstrap by heat, resampling entire heats; avoids treating overlapping windows as independent.
        rng=np.random.default_rng(SEED+len(slope_rows)+9000)
        hs=sorted(z.heat.unique()); boots=[]
        if len(hs)>=3:
            for _ in range(500):
                sampled=rng.choice(hs,size=len(hs),replace=True)
                pieces=[]
                for j,h in enumerate(sampled):
                    g=z[z.heat==h].copy(); g['boot_cluster']=f'{j}_{h}'; pieces.append(g)
                bd=pd.concat(pieces,ignore_index=True)
                dx=bd.flow_median_m3s-bd.groupby('boot_cluster').flow_median_m3s.transform('mean')
                dy=bd.best_lag_s-bd.groupby('boot_cluster').best_lag_s.transform('mean')
                den=float(np.dot(dx,dx))
                if den>0: boots.append(float(np.dot(dx,dy)/den))
        if boots:
            ci=np.quantile(boots,[.025,.5,.975])
            neg_share=float(np.mean(np.asarray(boots)<0))
        else: ci=[np.nan]*3; neg_share=np.nan
        # Leave-one-heat-out slope sign stability.
        loo=[]
        for h in hs:
            b=within_slope(z[z.heat!=h])
            loo.append(b)
        slope_rows.append({'scenario':scenario,'gas':gas,'n_informative_windows':len(z),'n_heats':len(hs),
            'within_heat_slope_s_per_m3s':beta,'cluster_boot_ci_low':ci[0],'cluster_boot_median':ci[1],
            'cluster_boot_ci_high':ci[2],'bootstrap_share_beta_negative':neg_share,
            'loo_negative_share':float(np.mean(np.asarray(loo)<0)) if loo else np.nan,
            'median_per_heat_spearman':per_df.rho.median() if len(per_df) else np.nan,
            'heats_with_negative_rho':int((per_df.rho<0).sum()) if len(per_df) else 0,
            'heats_with_positive_rho':int((per_df.rho>0).sum()) if len(per_df) else 0})
        for _,r in per_df.iterrows(): slope_boot_rows.append({'scenario':scenario,'gas':gas,**r.to_dict()})
slope_summary=pd.DataFrame(slope_rows)
slope_summary.to_csv(OUT/'tables/within_heat_flow_delay_slope.csv',index=False)
pd.DataFrame(slope_boot_rows).to_csv(OUT/'tables/per_heat_flow_delay_slopes.csv',index=False)

# Leave-one-heat-out prediction of window lag from flow: fixed mean vs linear flow vs monotone non-increasing isotonic.
model_rows=[]
for scenario in ['H0','H1']:
    for gas in GASES:
        z=windows[(windows.scenario==scenario)&(windows.gas==gas)&(windows.informative)&(windows.in_reference_flow_support)].copy()
        if scenario=='H0': z=z[z.heat!='6C02057']
        heats=sorted(z.heat.unique())
        for testh in heats:
            tr=z[z.heat!=testh]; te=z[z.heat==testh]
            if len(te)<1 or tr.heat.nunique()<3: continue
            ytr=tr.best_lag_s.to_numpy(); xtr=tr.flow_median_m3s.to_numpy().reshape(-1,1)
            yte=te.best_lag_s.to_numpy(); xte=te.flow_median_m3s.to_numpy().reshape(-1,1)
            pred_fixed=np.repeat(np.median(ytr),len(te))
            lm=LinearRegression().fit(xtr,ytr); pred_lin=lm.predict(xte)
            # Isotonic non-increasing (the problem's directional prior), trained only on other heats.
            iso=IsotonicRegression(increasing=False,out_of_bounds='clip').fit(xtr.ravel(),ytr)
            pred_iso=iso.predict(xte.ravel())
            for m,pred in [('fixed_median',pred_fixed),('free_linear',pred_lin),('nonincreasing_isotonic',pred_iso)]:
                model_rows.append({'scenario':scenario,'gas':gas,'test_heat':testh,'model':m,
                    'n_test_windows':len(te),'test_MAE_s':float(np.mean(np.abs(pred-yte))),
                    'test_RMSE_s':float(np.sqrt(np.mean((pred-yte)**2)))})
model_df=pd.DataFrame(model_rows)
model_df.to_csv(OUT/'tables/flow_delay_model_loho.csv',index=False)

# Spectral brightness proxy: 32 equal-width bands on 345-1015 nm, log of mean band brightness.
spec_rows=[]
if SPEC_FILE.exists():
    zsp=np.load(SPEC_FILE)
    grid=zsp['grid'].astype(float)
    mask=(grid>=345)&(grid<=1015)
    edges=np.linspace(345,1015,33)
    band_idx=[]
    for a,b in zip(edges[:-1],edges[1:]): band_idx.append(np.where(mask&(grid>=a)&(grid<b))[0])
    for scenario in ['H0','H1']:
        for h in HEATS:
            X=zsp[h].astype(np.float32)
            bm=np.column_stack([np.mean(X[:,ix],axis=1) for ix in band_idx])
            brightness=np.log(np.maximum(np.mean(bm,axis=1),1e-6))
            dt0=1.0 if (scenario=='H1' and h=='6C02057') else BASE_DT
            stride=max(1,int(round(1.0/dt0)))
            brightness=brightness[::stride]
            for gas in GASES:
                g=curves[h][f'{gas}_flow'].to_numpy(float)[::stride]
                dt=dt0*stride
                for sw in [5.0,9.0,15.0]:
                    x=smooth_derivative(brightness,dt,sw)
                    y=smooth_derivative(g,dt,sw)
                    trim=min(int(round(20/dt)),max(0,(len(x)-100)//4))
                    xx=x[trim:-trim] if trim else x; yy=y[trim:-trim] if trim else y
                    ll,rr=corr_lags(xx,yy,dt,MAX_LAG_S)
                    sm=peak_summary(ll,rr)
                    spec_rows.append({'scenario':scenario,'heat':h,'gas':gas,'smooth_s':sw,
                        'dt_s':dt,'n_points':len(xx),**sm})
spec_df=pd.DataFrame(spec_rows)
if len(spec_df): spec_df.to_csv(OUT/'tables/spectral_brightness_lag_diagnostics.csv',index=False)

# Reconstruct old/previous baseline summaries for comparison (not re-used as ground truth).
legacy_path=ROOT/'Q2_Figures_Package/07_Supporting_Data/q2_flow_lag_bins.csv'
if legacy_path.exists():
    legacy=pd.read_csv(legacy_path)
    legacy.to_csv(OUT/'tables/previous_60s_window_flow_bins_reference.csv',index=False)

# Save an explicit data-quality inventory.
quality=[]
for h,d in curves.items():
    for gas in GASES:
        a=d[f'{gas}_flow'].to_numpy(float)
        quality.append({'heat':h,'gas':gas,'n_rows':len(a),'zero_share':float(np.mean(a==0)),
            'negative_count':int(np.sum(a<0)),'min':float(np.min(a)),'median':float(np.median(a)),
            'max':float(np.max(a)),'sd':float(np.std(a)),'unique_values':int(pd.Series(a).nunique())})
pd.DataFrame(quality).to_csv(OUT/'tables/gas_signal_quality.csv',index=False)

# Figures.
# 1 full-series lag curves for all heats (H0), one gas per page.
for gas in GASES:
    fig,axes=plt.subplots(2,5,figsize=(18,7),sharex=True,constrained_layout=True)
    for ax,h in zip(axes.flat,HEATS):
        s=series_cache[('H0',h,gas)]
        ax.plot(s['lags'],s['corr'],linewidth=1.2)
        ax.axhline(0,linewidth=.7)
        ax.set_title(h,fontsize=9); ax.grid(True,alpha=.25)
        ax.set_xlabel('Gas response lag (s)'); ax.set_ylabel('Pearson r')
    fig.suptitle(f'Oxygen-flow change vs {gas} change — H0 (0.5 s/row)',fontsize=14)
    fig.savefig(OUT/f'figures/H0_full_lag_profiles_{gas}.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# 2 per-heat lag with bootstrap intervals, H0 and H1.
for scenario in ['H0','H1']:
    fig,axes=plt.subplots(1,2,figsize=(15,6),sharey=True,constrained_layout=True)
    for ax,gas in zip(axes,GASES):
        sub=full[(full.scenario==scenario)&(full.gas==gas)].set_index('heat').reindex(HEATS).reset_index()
        y=np.arange(len(sub))
        lo=sub.best_abs_lag_s-sub.bootstrap_lag_ci_low_s
        hi=sub.bootstrap_lag_ci_high_s-sub.best_abs_lag_s
        ax.errorbar(sub.best_abs_lag_s,y,xerr=[lo.clip(lower=0),hi.clip(lower=0)],fmt='o',capsize=3)
        ax.set_yticks(y,labels=sub.heat); ax.invert_yaxis(); ax.set_xlim(0,MAX_LAG_S)
        ax.set_xlabel('Best absolute-correlation lag (s)'); ax.set_title(gas); ax.grid(axis='x',alpha=.25)
    fig.suptitle(f'Per-heat exploratory delay and paired-block bootstrap interval — {scenario}',fontsize=14)
    fig.savefig(OUT/f'figures/{scenario}_per_heat_delay_bootstrap.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# 3 flow-delay windows scatter.
for scenario in ['H0','H1']:
    fig,axes=plt.subplots(1,2,figsize=(14,5.5),constrained_layout=True)
    for ax,gas in zip(axes,GASES):
        sub=windows[(windows.scenario==scenario)&(windows.gas==gas)]
        for h,g in sub.groupby('heat'):
            ax.scatter(g.flow_median_m3s,g.best_lag_s,s=24,alpha=.65,label=h)
        ax.axvspan(flow_lo,flow_hi,alpha=.08)
        ax.set_xlabel('Median oxygen flow in 120-s window (m³/s)')
        ax.set_ylabel('Exploratory local lag (s)'); ax.set_title(gas); ax.grid(True,alpha=.25)
        ax.legend(fontsize=7,ncol=2,frameon=False)
    fig.suptitle(f'Flow vs local lag (120-s windows, 60-s stride) — {scenario}',fontsize=14)
    fig.savefig(OUT/f'figures/{scenario}_flow_vs_local_lag.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# 4 flow-bin medians informative windows.
for scenario in ['H0','H1']:
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for ax,gas in zip(axes,GASES):
        z=bins[(bins.scenario==scenario)&(bins.gas==gas)&(bins.window_set=='informative')]
        x=np.arange(3); med=z.set_index('flow_bin').reindex(['low','middle','high']).lag_median_s.to_numpy()
        q1=z.set_index('flow_bin').reindex(['low','middle','high']).lag_q25_s.to_numpy()
        q3=z.set_index('flow_bin').reindex(['low','middle','high']).lag_q75_s.to_numpy()
        ax.errorbar(x,med,yerr=[med-q1,q3-med],fmt='o-',capsize=4)
        ax.set_xticks(x,labels=['Low flow','Middle flow','High flow'])
        ax.set_ylabel('Lag (s)'); ax.set_title(gas); ax.grid(True,alpha=.25)
    fig.suptitle(f'Flow-bin lag medians and IQR — {scenario}',fontsize=14)
    fig.savefig(OUT/f'figures/{scenario}_flow_bin_lag.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# 5 H0/H1 comparison for fifth furnace.
fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
for ax,gas in zip(axes,GASES):
    sub=full[(full.heat=='6C02057')&(full.gas==gas)].set_index('scenario').reindex(['H0','H1'])
    ax.bar(['H0 0.5 s','H1 1.0 s'],sub.best_abs_lag_s.to_numpy())
    ax.set_ylabel('Best absolute-correlation lag (s)'); ax.set_title(gas); ax.grid(axis='y',alpha=.25)
fig.suptitle('Sixth? No — Heat 6C02057: time-scale sensitivity of estimated lag',fontsize=13)
fig.savefig(OUT/'figures/6C02057_H0_H1_lag_sensitivity.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# 6 spectral brightness delay profiles summary.
if len(spec_df):
    fig,axes=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
    for ax,gas in zip(axes,GASES):
        sub=spec_df[(spec_df.scenario=='H0')&(spec_df.gas==gas)]
        for sw,g in sub.groupby('smooth_s'):
            ax.scatter(g.best_abs_lag_s,g.best_abs_r,label=f'{sw:g}s smooth',alpha=.7)
        ax.set_xlabel('Best absolute-correlation lag (s)'); ax.set_ylabel('Peak |r|')
        ax.set_title(f'Spectral brightness → {gas}'); ax.grid(True,alpha=.25); ax.legend(frameon=False)
    fig.suptitle('Supplementary optical-to-gas lag diagnostic (H0)',fontsize=14)
    fig.savefig(OUT/'figures/H0_spectral_to_gas_lag_diagnostic.png',dpi=200,bbox_inches='tight'); plt.close(fig)

# Summary stats for report.
summary={}
for scenario in ['H0','H1']:
    summary[scenario]={}
    for gas in GASES:
        z=full[(full.scenario==scenario)&(full.gas==gas)]
        w=windows[(windows.scenario==scenario)&(windows.gas==gas)]
        wi=w[w.informative]
        summary[scenario][gas]={
            'median_full_lag_s':float(z.best_abs_lag_s.median()),
            'median_peak_abs_r':float(z.best_abs_r.abs().median()),
            'heats_with_BH_q_le_0_10':int((z.BH_q<=.10).sum()),
            'median_bootstrap_ci_width_s':float((z.bootstrap_lag_ci_high_s-z.bootstrap_lag_ci_low_s).median()),
            'informative_windows':int(len(wi)),'informative_heats':int(wi.heat.nunique()),
            'informative_lag_median_s':float(wi.best_lag_s.median()) if len(wi) else None,
            'informative_flow_median_m3s':float(wi.flow_median_m3s.median()) if len(wi) else None}

# Report flow slopes.
summary['flow_delay_slopes']=slope_summary.replace({np.nan:None}).to_dict(orient='records')
summary['common_flow_support']={'q01_m3s':flow_lo,'q99_m3s':flow_hi,'tertile_cut_1_m3s':float(qcuts[0]),'tertile_cut_2_m3s':float(qcuts[1])}
summary['spectral']={}
if len(spec_df):
    for gas in GASES:
        z=spec_df[(spec_df.scenario=='H0')&(spec_df.gas==gas)]
        summary['spectral'][gas]={'median_peak_abs_r':float(z.best_abs_r.abs().median()),
            'median_lag_s':float(z.best_abs_lag_s.median()),
            'boundary_peak_count':int(z.boundary.sum()),'n_heat_smoothing_pairs':int(len(z))}
(OUT/'tables/summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

# Human-readable report.
lines=[]
lines += ['# APMCM 2026 B题第二问：延迟识别专项实验', '',
'版本：Q2-delay-v1.0；数据：十炉连续曲线与 harmonized_spectra.npz；实验日期：2026-11-26。', '',
'## 1. 本轮实验目标',
'检验供氧变化到 CO₂、CO 烟气流量变化之间是否存在可重复的正向时间延迟；进一步检验延迟是否随供氧流量增加而缩短。对第五炉 6C02057 保留双时间尺度情景。', '',
'## 2. 时间尺度情景',
'- H0（主情景）：全部记录按 0.5 s/行。',
'- H1（敏感性）：仅 6C02057 按 1.0 s/行，其余仍为 0.5 s/行。',
'- H1 不是数据修复结论；当前附件不能识别第五炉的绝对采样时间尺度。任何涉及秒、m³/s、积分和延迟的结果均需并列报告。', '',
'## 3. 方法',
'- 上游参照：由累计耗氧量在 10 s 居中窗口内做局部线性斜率，得到非负氧枪流量估计；再对氧流与气体流量分别做 5 s Savitzky–Golay 平滑并求一阶变化率。',
'- 延迟定义：对氧流变化率与气体变化率计算非负滞后互相关；正滞后表示气体信号较晚。搜索范围 0–60 s。报告最大绝对相关、正相关峰、负相关峰、零滞后相关和峰宽。',
'- 不确定性：对配对序列做 120 s 移动块重抽样，给出峰值滞后的 95%百分位区间；用循环移位负对照保持单条序列的自相关结构，并对每个情景的 20 个炉次×气体检验做 BH 多重检验校正。',
'- 流量条件分析：使用 120 s 窗口、60 s 步长、0–40 s 局部滞后搜索。只有流量标准差≥0.15 m³/s、峰值 |r|≥0.25 且峰不在边界的窗口标记为“信息性窗口”。此阈值是探索性筛选，不是物理阈值。',
'- 流量分箱边界按 H0 下除第五炉外的九炉窗口流量分布三等分确定；同一边界用于 H0/H1，避免因情景切换重定义分箱。第五炉 H0 若超出九炉共同支持域，单独标记，不用于通用斜率估计。',
'- 辅助光谱诊断：在 345–1015 nm 取 32 个等宽波段均值，构造对数平均亮度，按 1 Hz 取样，在 5/9/15 s 平滑尺度下计算光谱亮度变化到烟气变化的互相关。该项仅作辅助诊断，不替代独立光谱状态建模。', '',
'## 4. 核心结果',
'### 4.1 全炉互相关',
'逐炉结果见 `tables/full_series_lag_diagnostics.csv`。互相关峰值是探索性有效时差，不是已知物理运输时间。若峰值相关弱、峰宽大、峰落在边界或循环移位检验不支持，则该炉该气体不能给出稳定延迟点估计。', '',
'### 4.2 流量—延迟关系',
'逐窗口数据见 `tables/local_window_lag_diagnostics.csv`；分箱统计见 `tables/flow_bin_lag_summary.csv`；炉级聚类 bootstrap 斜率见 `tables/within_heat_flow_delay_slope.csv`。斜率为负才符合“流量越大、延迟越短”的方向。只有在不同炉次方向稳定、区间支持域内复现且区间排除零时，才可称为支持该规律。', '',
'### 4.3 第五炉时间尺度敏感性',
'见 `figures/6C02057_H0_H1_lag_sensitivity.png` 及逐炉表。H1 会同时改变第五炉的物理时间轴、氧流量尺度、局部窗口数量和延迟秒数；它不是简单把最终 lag 乘 2。H0 作为主结果，H1 只作敏感性。', '',
'### 4.4 负对照与区间',
'循环移位 p 值和 BH q 值在逐炉表中。移动块 bootstrap 区间若宽、峰值位置不稳定，说明单炉信号不能精确识别时滞。循环移位仅为保持自相关的近似负对照；非平稳序列下不构成严格随机化检验。', '',
'### 4.5 光谱到烟气辅助诊断',
'见 `tables/spectral_brightness_lag_diagnostics.csv`。该特征是光谱亮度代理，不等于第三问最终冻结的光谱状态；不能把供氧到烟气时差替代光谱到烟气时差。', '',
'## 5. 模型升级判断',
'本轮不预设非递增时滞模型胜出。若流量—延迟斜率在炉级 bootstrap 中不稳定，或 LOHO 中流量线性/非递增模型不能优于固定中位数基准，应保留固定有效时差或报告弱识别，而不是强制施加单调约束。', '',
'## 6. 重要限制',
'- 氧枪瞬时流量由累计耗氧反推，无独立流量传感器真值；其误差会传递到延迟估计。',
'- 烟气信号按题面数值口径作为瞬时体积流量 m³/s；如果官方/厂方定义不同，所有流量和积分解释需调整。',
'- 互相关峰表示两个变化序列的统计对齐，不可拆分炉内反应、烟道输运、混合和仪器响应。',
'- 十炉样本很少，窗口重叠且同炉自相关强；所有炉级统计均需谨慎，不把窗口数当独立炉次样本数。',
'- 若 CO₂、CO 的峰值方向在炉次间翻转，不能用单一正向滞后解释全部过程。',
'- 第五炉绝对时间尺度不可识别，H0/H1 需并列保留。', '',
'## 7. 文件索引',
'- `tables/full_series_lag_diagnostics.csv`：全炉、全气体、H0/H1 的峰值、bootstrap 区间、循环移位 p 与 BH q。',
'- `tables/local_window_lag_diagnostics.csv`：逐窗口流量、相关峰、延迟和信息性标记。',
'- `tables/flow_bin_lag_summary.csv`：流量分箱统计。',
'- `tables/within_heat_flow_delay_slope.csv`：炉级聚类 bootstrap 的流量—延迟斜率。',
'- `tables/per_heat_flow_delay_slopes.csv`：逐炉流量—延迟相关。',
'- `tables/flow_delay_model_loho.csv`：固定、自由线性、非递增等模型的留一炉预测误差。',
'- `tables/spectral_brightness_lag_diagnostics.csv`：辅助光谱亮度—烟气延迟。',
'- `tables/gas_signal_quality.csv`：烟气序列基本质量统计。',
'- `code/run_q2_delay_experiment.py`：本轮完整复现脚本。', '',
'## 8. 当前结论状态',
'本实验用于判断延迟关系是否可识别，不保证得出题面预期的单调规律。最终报告应以数据支持为准；若证据不够，应明确写“现有十炉数据不足以稳定识别流量依赖的延迟”，并给出逐炉探索性时差范围和失败原因。']
(OUT/'README_Q2_Delay_Experiment.md').write_text('\n'.join(lines),encoding='utf-8')

# Copy input curves and provenance.
for h in HEATS:
    shutil.copy2(CURVE_DIR/f'curves_{h}.csv', OUT/'inputs'/f'curves_{h}.csv')
shutil.copy2(SPEC_FILE, OUT/'inputs/harmonized_spectra.npz')

# Copy this script into package.
shutil.copy2(Path(__file__), OUT/'code/run_q2_delay_experiment.py')

# Write source fingerprints.
fp=[]
for p in sorted((OUT/'inputs').iterdir()):
    if p.is_file():
        sha=hashlib.sha256(p.read_bytes()).hexdigest()
        fp.append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha})
pd.DataFrame(fp).to_csv(OUT/'tables/input_sha256.csv',index=False)

# Build ZIP.
zip_path=ROOT.parent/'Q2_Delay_Identification_Experiment_rerun.zip'
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
    for p in sorted(OUT.rglob('*')):
        if p.is_file(): zf.write(p,arcname=str(p.relative_to(OUT)))

print('OUTPUT',zip_path,zip_path.stat().st_size)
print('\nFULL SUMMARY')
print(full[['scenario','heat','gas','best_abs_lag_s','best_abs_r','bootstrap_lag_ci_low_s','bootstrap_lag_ci_high_s','circular_shift_p','BH_q','boundary']].round(3).to_string(index=False))
print('\nSLOPE SUMMARY')
print(slope_summary.round(3).to_string(index=False))
print('\nFLOW BINS informative')
print(bins[bins.window_set=='informative'].round(2).to_string(index=False))
print('\nLOHO MODEL SUMMARY')
if len(model_df): print(model_df.groupby(['scenario','gas','model']).test_MAE_s.mean().round(2).to_string())
print('\nSPECTRAL SUMMARY H0')
if len(spec_df): print(spec_df[spec_df.scenario=='H0'].groupby(['gas','smooth_s']).agg(median_abs_r=('best_abs_r',lambda x: np.median(np.abs(x))),median_lag=('best_abs_lag_s','median'),boundary_count=('boundary','sum')).round(3).to_string())
print('\nSUMMARY_JSON')
print(json.dumps(summary,indent=2,ensure_ascii=False))
