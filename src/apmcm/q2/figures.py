import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from .reconstruct import ROOT,TAB

OUT=ROOT/'figures'
TEAL='#007F86'
ORANGE='#BC5B28'
BLUE='#35618C'
DARK='#24394A'
plt.rcParams.update({'figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,
                     'axes.labelcolor':DARK,'text.color':DARK,'font.size':9,'axes.titleweight':'bold',
                     'axes.titlelocation':'left','savefig.dpi':165,'axes.axisbelow':True})


def finish(fig,name):
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight',pad_inches=.14)
    plt.close(fig)


def figures():
    OUT.mkdir(exist_ok=True)
    curves=pd.read_csv(TAB/'curves.csv.gz')
    for heat,g in curves[curves.scenario=='H0'].groupby('heat'):
        fig,axs=plt.subplots(3,2,figsize=(10,6.9),layout='constrained',sharex=True)
        a=int(g.attachment.iloc[0])
        fig.suptitle(f'Attachment {a:02d}    /    Heat {heat}',x=.03,ha='left',fontsize=15,fontweight='bold')
        for ax,var,title,unit in zip(axs.flat,['O','q','CO2','CO','E_CO2','E_CO'],
            ['Cumulative oxygen','Lance flow','CO₂ flow','CO flow','Cumulative CO₂','Cumulative CO'],
            ['m³','m³/s','m³/s','m³/s','m³','m³']):
            if var+'_raw' in g:
                ax.scatter(g.t,g[var+'_raw'],s=4,color='#ACB5BF',rasterized=True,label='Recorded')
            ax.plot(g.t,g[var],lw=1.3,color=TEAL if var not in ['CO','E_CO'] else ORANGE,label='Reconstruction')
            ax.set(title=title,ylabel=unit)
            ax.grid(axis='y',alpha=.2)
        for ax in axs[-1]: ax.set_xlabel('Time (s)')
        axs[0,0].legend(frameon=False,fontsize=8)
        finish(fig,f'heat_{a:02d}_{heat}')
    folds=pd.read_csv(TAB/'delay_folds.csv').query("scenario=='H0' and window==8")
    fig,axes=plt.subplots(1,2,figsize=(9,3.25),layout='constrained')
    for ax,(gas,g) in zip(axes,folds.groupby('gas',sort=True)):
        p=g.pivot(index='heat',columns='model',values='mae')
        v=(p.decreasing-p.fixed).sort_index()
        ax.barh(range(10),v,color=[TEAL if n<0 else ORANGE for n in v])
        ax.axvline(0,c=DARK,lw=.7)
        ax.set(yticks=range(10),yticklabels=v.index,title=gas.replace('2','₂'),xlabel='MAE difference (m³/s)')
        ax.grid(axis='x',alpha=.18)
        ax.invert_yaxis()
    finish(fig,'delay_paired')
    glob=pd.read_csv(TAB/'global_profiles.csv').query("scenario=='H0'")
    fig,axs=plt.subplots(1,2,figsize=(9,3.1),layout='constrained',sharey=True)
    for ax,(gas,g) in zip(axs,glob.groupby('gas')):
        for heat,v in g.groupby('heat'):
            ax.plot(v.lag_s,v.r,color=BLUE,alpha=.28,lw=1)
        med=g.groupby('lag_s').r.median()
        ax.plot(med.index,med.values,color=DARK,lw=2,label='Median by lag')
        ax.axhline(0,c='#95A2AD',lw=.7)
        ax.set(title=gas.replace('2','₂'),xlabel='Gas lag (s)',ylabel='Signed correlation')
        ax.legend(frameon=False,fontsize=8)
    finish(fig,'global_profiles')
    bins=pd.read_csv(TAB/'empirical_flow_bins.csv')
    fig,axes=plt.subplots(1,2,figsize=(9,2.9),layout='constrained',sharey=True)
    for ax,(scen,g) in zip(axes,bins.groupby('scenario')):
        for j,(gas,sub) in enumerate(g.groupby('gas')):
            x=np.arange(3)+(j-.5)*.14
            ax.errorbar(x,sub.lag_median,yerr=[sub.lag_median-sub.lag_q25,sub.lag_q75-sub.lag_median],
                        fmt='o-',capsize=3,lw=1.5,color=[ORANGE,TEAL][j],label=gas.replace('2','₂'))
        ax.set(title=scen,xticks=range(3),xticklabels=['Low','Middle','High'],ylabel='Peak lag (s)',xlabel='Oxygen flow tertile',ylim=(-1,31))
        ax.grid(axis='y',alpha=.2)
        ax.legend(frameon=False,ncol=2)
    finish(fig,'local_flow_bins')
    fig,axes=plt.subplots(1,2,figsize=(9,3),layout='constrained')
    for scen,h in curves[curves.attachment==5].groupby('scenario'):
        for ax,var,unit in zip(axes,['q','E_CO'],['Lance flow (m³/s)','Cumulative CO (m³)']):
            ax.plot(h.t,h[var],color=TEAL if scen=='H0' else ORANGE,lw=1.2,label=scen)
            ax.set(xlabel='Time (s)',title=unit)
            ax.legend(frameon=False)
    finish(fig,'fifth_timebase')
    a=pd.read_csv(TAB/'associations.csv').query("scenario=='H0'")
    fig,ax=plt.subplots(figsize=(8,2.7),layout='constrained')
    for j,var in enumerate(['CO2','CO','E_CO2','E_CO']):
        v=a[a.outcome==var]
        for col,dx,c in [('rho_level',-.13,BLUE),('rho_blocks',.13,ORANGE)]:
            ax.scatter(j+dx+np.linspace(-.05,.05,len(v)),v[col],s=13,alpha=.55,c=c)
            ax.plot([j+dx-.10,j+dx+.10],[v[col].median()]*2,c=c,lw=3)
    ax.set(xticks=range(4),xticklabels=['CO₂ flow','CO flow','Cumulative CO₂','Cumulative CO'],ylabel='Spearman correlation',ylim=(-1.06,1.12))
    ax.axhline(0,c=DARK,lw=.6)
    ax.legend(handles=[Line2D([0],[0],color=BLUE,label='Full levels'),Line2D([0],[0],color=ORANGE,label='30 s blocks')],loc='lower right',frameon=False)
    finish(fig,'oxygen_associations')


def dose_response():
    all_data=pd.read_csv(TAB/'curves.csv.gz')
    results=[]
    fig,axes=plt.subplots(2,2,figsize=(8.5,4.2),layout='constrained',sharex=True)
    levels=np.arange(0,3401,50.)
    for ax,var,label in zip(axes.flat,['CO2','CO','E_CO2','E_CO'],['CO₂ flow (m³/s)','CO flow (m³/s)','Cumulative CO₂ (m³)','Cumulative CO (m³)']):
        matrix=[]
        for heat,g in all_data[all_data.scenario=='H0'].groupby('heat'):
            x=(g.O-g.O.iloc[0]).to_numpy()
            y=g[var].to_numpy()
            ax.plot(x,y,c=BLUE,lw=.6,alpha=.25)
            v=np.interp(levels,x,y,left=np.nan,right=np.nan)
            matrix.append(v)
            for o,z in zip(levels,v):
                if np.isfinite(z):
                    results.append(dict(heat=heat,outcome=var,oxygen_net=o,value=z))
        ax.plot(levels,np.nanmedian(matrix,axis=0),c=ORANGE,lw=1.8,label='Median across available heats')
        ax.set(title=label)
        ax.grid(axis='y',alpha=.18)
    for ax in axes[-1]: ax.set_xlabel('Net cumulative oxygen (m³)')
    finish(fig,'dose_response')
    pd.DataFrame(results).to_csv(TAB/'dose_response.csv',index=False)


if __name__=="__main__":
    figures()
    dose_response()
