"""Scoped Q1 diagnostics. Run from the project root: python -m src.q1_peak_study."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from python_calamine import CalamineWorkbook
from scipy.signal import find_peaks
from scipy.ndimage import maximum_filter1d, label
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from .common import CURRENT_DATA, RAW

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/tables"
GRID = np.linspace(345,1015,1200)
CAP = 65535
SEED = 17291

def load_raw(row):
    path = RAW/row.file
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row.sha256
    rows = CalamineWorkbook.from_path(str(path)).get_sheet_by_index(0).iter_rows()
    w = np.asarray(next(rows)[1:-3],float)
    x = np.asarray([r[1:-3] for r in rows],dtype=np.float32)
    return w,x

def peaks(y,w):
    # Half-prominence width in nm, finite-span prominence, with a local brightness floor.
    pos,props = find_peaks(y,prominence=250,width=(None,None),wlen=31)
    widths=np.interp(props["right_ips"],np.arange(len(w)),w)-np.interp(props["left_ips"],np.arange(len(w)),w)
    return pos[(widths<=1.5)&(props["prominences"]>=np.maximum(250,.065*y[pos]))&(y[pos]<CAP)]

def negative_peaks(y,w):
    pos,props=find_peaks(-y,prominence=250,width=(None,None),wlen=31)
    axis=np.arange(len(w))
    widths=np.interp(props["right_ips"],axis,w)-np.interp(props["left_ips"],axis,w)
    return pos[(widths<=1.5)&(props["prominences"]>=np.maximum(250,.065*y[pos]))]

def shape_mask(x,w):
    out=np.zeros(x.shape,bool)
    for i,y in enumerate(x):
        out[i,peaks(y,w)]=True
    return out

def supported(mask):
    near=maximum_filter1d(mask.astype('uint8'),3,axis=1).astype(bool)
    both=np.zeros_like(mask)
    both[1:-1]=near[:-2]&near[2:]
    either=np.zeros_like(mask)
    either[1:-1]=near[:-2]|near[2:]
    return both,either

def isolated(y,i,centers,either):
    j=centers[(centers>=4)&(centers<y.shape[1]-4)]
    if not len(j): return j
    neighbors=(y[i-1]+y[i+1])/2
    # Remove local common-mode brightness changes using flanks outside the candidate.
    shifts=np.array([-4,-3,3,4])
    curr=y[i,j[:,None]+shifts]
    prev=neighbors[j[:,None]+shifts]
    scale=np.median(curr/np.maximum(prev,1),axis=1)
    ry=y[i,j]-scale*neighbors[j]
    noise=np.maximum(1,1.4826*np.median(abs(curr-scale[:,None]*prev),axis=1))
    return j[(~either[i,j])&(ry>np.maximum(6*noise,.08*y[i,j]))]

def support_stats(x):
    cap=x==CAP
    counts=cap.sum(1)
    left=np.pad(cap[:,:-1],((0,0),(1,0)))
    right=np.pad(cap[:,1:],((0,0),(0,1)))
    blocks,n=label(cap)
    return dict(cap=int(cap.sum()),cap_frames=int((counts>0).sum()),cap_adjacent=int((cap&(left|right)).sum()),
                largest_cap_patch=int(np.bincount(blocks.ravel())[1:].max()) if n else 0)

def write(name,rows):
    df=pd.DataFrame(rows)
    df.to_csv(OUT/("q1_peak_"+name+".csv"),index=False)
    return df

def injection(h,w,x,rng):
    # Balanced positive/untouched controls. Copies are modified for one trial only.
    outputs=[]
    for level in [.12,.4]:
        for n in range(160):
            i=rng.integers(2,len(x)-2)
            j=rng.integers(10,x.shape[1]-10)
            patch=x[i-1:i+2,j-20:j+21].copy() if 20<=j<x.shape[1]-20 else None
            if patch is None or patch.max()>50000: continue
            bands=w[j-20:j+21]
            add=level*max(2500.,float(patch[1,20]))
            for yes in [0,1]:
                p=patch.copy()
                if yes: p[1,20]+=add
                masks=shape_mask(p,bands); _,either=supported(masks)
                candidates=np.flatnonzero(masks[1]); joint=isolated(p,1,candidates,either)
                for model,flag in [('shape',masks[1,20]),('joint',20 in joint)]:
                    outputs.append(dict(heat=h,amplitude=level,model=model,injected=yes,flag=int(flag)))
    return outputs

def projection(train):
    # Equal heat weight and within-heat centering for the offline reconstruction basis.
    centered=np.concatenate([(v-v.mean(0))/np.sqrt(len(v)) for v in train])
    from sklearn.utils.extmath import randomized_svd
    _,_,vt=randomized_svd(centered,n_components=12,n_iter=4,random_state=SEED)
    return vt

def patch_scores(pred,true,w,bg):
    h0=np.maximum(true-bg,0);h1=np.maximum(pred-bg,0)
    area=np.trapezoid(h0,w)
    return dict(mae=np.mean(abs(pred-true)),mse=np.mean((pred-true)**2),
                height_err=abs(h1.max()-h0.max())/max(h0.max(),1),
                area_err=abs(np.trapezoid(h1,w)-area)/max(area,1),
                location_nm=abs(w[h1.argmax()]-w[h0.argmax()]))

def study():
    manifest=pd.read_csv(CURRENT_DATA/"manifest.csv")
    data={};sample={};aud=[];inj=[];scores=[];clips=[];examples={}
    for row in manifest.itertuples():
        h=row.heat;w,x=load_raw(row)
        sample[h]=np.array([np.interp(GRID,w,y) for y in x[::20]])
        data[h]=(w,x)
    for idx,(h,(w,x)) in enumerate(data.items()):
        rng=np.random.default_rng(SEED+idx)
        mask=shape_mask(x,w);both,either=supported(mask)
        flags=np.zeros(x.shape,bool)
        for i in range(1,len(x)-1):
            flags[i,isolated(x,i,np.flatnonzero(mask[i]),either)]=True
        trough=np.zeros_like(mask)
        for i,y in enumerate(x):
            trough[i,negative_peaks(y,w)]=True
        t_both,t_either=supported(trough)
        neg=0
        for i in range(1,len(x)-1):
            j=np.flatnonzero(trough[i]&~t_either[i]);j=j[(j>=4)&(j<x.shape[1]-4)]
            if not len(j):continue
            # A negative residual must be large compared with the same shoulder noise floor.
            nei=(x[i-1]+x[i+1])/2
            jj=j[:,None]+np.array([-4,-3,3,4]);yg=x[i,jj];ng=nei[jj]
            sc=np.median(yg/np.maximum(ng,1),axis=1)
            noise=np.maximum(1,1.4826*np.median(abs(yg-sc[:,None]*ng),axis=1))
            neg+=np.count_nonzero(sc*nei[j]-x[i,j]>np.maximum(6*noise,.08*x[i,j]))
        cap=x==CAP
        aud.append(dict(heat=h,frames=len(x),cells=x.size,eligible_centers=int(mask[1:-1].sum()),
                        recurring=int((both&mask).sum()),no_adjacent=int((mask[1:-1]&~either[1:-1]).sum()),
                        joint_flags=int(flags.sum()),negative_candidates=int(trough[1:-1].sum()),negative_both=int((t_both&trough).sum()),negative_flags=neg,**support_stats(x)))
        inj.extend(injection(h,w,x,rng))
        vt=projection([v for g,v in sample.items() if g!=h])
        # On native channels evaluate common basis, then fit scores on visible channels only.
        common=(w>=GRID[0])&(w<=GRID[-1])
        wn=w[common];basis=np.column_stack([np.interp(wn,GRID,b) for b in vt])
        off=np.flatnonzero(common)[0]
        # Static mean from other times in this heat; exclude all target frames +/- 1.
        times=np.arange(5,len(x)-5,17)
        allowed=np.ones(len(x),bool)
        for dt in [-1,0,1]: allowed[times+dt]=False
        avg=x[allowed][:,common].mean(0).astype(float)
        own=sample[h][allowed[::20]]
        local_v=projection([own])
        opposite=projection([v for g,v in sample.items() if g[:2]!=h[:2]])
        extras={key:np.column_stack([np.interp(wn,GRID,b) for b in v]) for key,v in
                [('within_rank12',local_v),('opposite_grid_rank12',opposite)]}
        for i in times:
            # Select test region from the two neighboring spectra before opening target centers.
            c=peaks((x[i-1]+x[i+1])/2,w)
            c=c[(w[c]>350)&(w[c]<1010)&both[i,c]]
            if not len(c):continue
            for j in rng.choice(c,min(3,len(c)),replace=False):
                js=j+np.arange(-6,7);inside=np.arange(3,10);ex=np.array([0,1,2,10,11,12])
                if x[i-1:i+2,js].max()>=CAP:continue
                ww=w[js];target=x[i,js].astype(float)
                prior=(x[i-1,js]+x[i+1,js]).astype(float)/2
                # Fit the local brightness factor from visible shoulders.
                gain=np.median(target[ex]/np.maximum(prior[ex],1))
                local=prior*gain
                background=np.interp(ww,[ww[2],ww[10]],[target[2],target[10]])
                core=target[inside];coords=ww[inside];bg=background[inside]
                if np.max(core-bg)<250:continue  # scoring stratum only
                known=np.ones(len(wn),bool);known[js[inside]-off]=False
                y=x[i,common]-avg
                pred={'wavelength_PCHIP':PchipInterpolator(ww[ex],target[ex])(coords),
                      'time_linear':prior[inside], 'shoulder_template':local[inside]}
                for r in [3,12]:
                    phi=basis[:,:r];coef=np.linalg.lstsq(phi[known],y[known],rcond=None)[0]
                    pred['shared_rank'+str(r)]=(avg+phi@coef)[js[inside]-off]
                for name,phi in extras.items():
                    coef=np.linalg.lstsq(phi[known],y[known],rcond=None)[0]
                    pred[name]=(avg+phi@coef)[js[inside]-off]
                shoulders=avg[js-off]
                pred['heat_mean_scaled']=np.median(target[ex]/np.maximum(shoulders[ex],1))*shoulders[inside]
                for alpha in [.25,.5,.75]:
                    pred['within_blend_'+str(alpha)]=(1-alpha)*pred['within_rank12']+alpha*pred['shoulder_template']
                    pred['blend_'+str(alpha)]=(1-alpha)*pred['shared_rank12']+alpha*pred['shoulder_template']
                # Local penalized least squares: adapt curvature on the supported ridge.
                d2=np.diff(np.eye(13),n=2,axis=0);obs=np.zeros(13);obs[ex]=1
                for name,weights in [('uniform',np.ones(11)),('adaptive',np.where(abs(np.arange(1,12)-6)<=3,.12,1.))]:
                    cov=np.diag(obs+.8)+4*d2.T@(weights[:,None]*d2)
                    z=solve(cov,obs*target+.8*local,assume_a='pos')
                    pred[name]=z[inside]
                for name,pp in pred.items():
                    scores.append(dict(heat=h,frame=i,center_nm=w[j],model=name,**patch_scores(pp,core,coords,bg)))
                # Artificial censoring at 70% of local peak prominence. The uncensored
                # neighbors and current unclipped shoulders stay visible in this test.
                ceiling=bg.mean()+.7*(core.max()-bg.mean()); hidden=core>ceiling
                if hidden.any():
                    for name,pp in [('clipped',np.minimum(core,ceiling)),
                                    ('template_lower_bound',np.where(hidden,np.maximum(ceiling,local[inside]),core))]:
                        clips.append(dict(heat=h,method=name,mae=np.mean(abs(pp[hidden]-core[hidden])),n=hidden.sum()))
                if h not in examples and max(core-bg)>1000:
                    examples[h]=(ww,target,prior,local,pred['shared_rank3'],pred['adaptive'],int(i),int(j))
    ad=write('audit',aud);d=write('mask_tests',scores);ic=write('injection',inj);cl=write('pseudoclip',clips)
    print(ad.to_string(index=False))
    # Freeze a single blend per heat on early patches, evaluate on disjoint late frames.
    split=[]
    for h,events in d.groupby('heat'):
        n=len(data[h][1])
        early=events[events.frame<n*.5-1]
        test=events[events.frame>n*.5+1]
        choices=['shared_rank12','within_rank12','shoulder_template']+[name+str(a) for name in ['blend_','within_blend_'] for a in [.25,.5,.75]]
        chosen=early[early.model.isin(choices)].groupby('model').mae.mean().idxmin()
        for model in ['shared_rank12','within_rank12','opposite_grid_rank12','shoulder_template','adaptive',chosen]:
            tmp=test[test.model==model]
            split.append(dict(heat=h,selected=chosen,model=model,n=len(tmp),mae=tmp.mae.mean(),mse=tmp.mse.mean(),
                              height_err=tmp.height_err.mean(),area_err=tmp.area_err.mean(),location_nm=tmp.location_nm.mean()))
    write('blocked',pd.DataFrame(split).drop_duplicates())
    print(d.groupby(['heat','model'])[['mae','mse','height_err','area_err','location_nm']].mean().groupby('model').mean().to_string())
    print('n masked',d.groupby('model').size().to_dict())
    print(ic.groupby(['amplitude','model','injected']).flag.agg(['sum','count']))
    print(cl.groupby(['heat','method']).mae.mean().groupby('method').mean())
    fig,axes=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    h=next(iter(examples));ww,y,prior,local,low,adap,i,j=examples[h]
    ax=axes[0];ax.plot(ww,y,'o-',ms=3,label="Observed");ax.plot(ww[3:10],low,'--',label="Rank 3");ax.plot(ww[3:10],adap,label="Adaptive")
    ax.axvspan(ww[3],ww[9],color='grey',alpha=.12);ax.legend(fontsize=7)
    ax.set(title=f"Masked peak | {h}, {i*.5:g} s",xlabel="Wavelength (nm)",ylabel="Recorded intensity")
    maxh=ad.loc[ad.cap.idxmax(),'heat'];w,x=data[maxh];i=int(np.argmax((x==CAP).sum(1)));ids=np.where(x[i]==CAP)[0]
    j=ids[len(ids)//2];js=np.arange(max(0,j-30),min(len(w),j+31));ax=axes[1]
    for shift in [-1,0,1]:ax.plot(w[js],x[i+shift,js],label=f"{.5*(i+shift):g} s")
    ax.axhline(CAP,ls=':',color='k');ax.legend(fontsize=7);ax.set(title=f"Upper-code plateau | {maxh}",xlabel="Wavelength (nm)")
    axes[2].scatter(ad.heat,100*ad.recurring/ad.eligible_centers,s=30)
    axes[2].tick_params(axis='x',labelrotation=75,labelsize=7)
    axes[2].set(title="Two-sided peak recurrence",ylabel="Share of narrow candidates (%)",xlabel="Heat")
    fig.savefig(ROOT/"outputs/figures/Q1_peak_evidence.png",dpi=165);plt.close(fig)

if __name__=='__main__':
    with threadpool_limits(1):
        study()
