import os,glob
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/mnt/data/APMCM_Full_Experiments';C=R+'/cleaned_scenarios';figdir=R+'/figures';os.makedirs(figdir,exist_ok=True)
for f in sorted(glob.glob(C+'/*_scenario.npz')):
 z=np.load(f);h=os.path.basename(f).split('_')[0];A=z['scenario'];wl=z['wl'];t=(z['seq']-1)*.5
 ti=np.arange(0,len(t),10);wi=np.arange(0,len(wl),10);T,W=np.meshgrid(t[ti],wl[wi],indexing='ij');Z=A[np.ix_(ti,wi)]
 fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d');ax.plot_surface(T,W,Z,cmap='viridis',linewidth=0,antialiased=False,rstride=1,cstride=1)
 ax.set_xlabel('Time (s)');ax.set_ylabel('Wavelength (nm)');ax.set_zlabel('Intensity (raw units)');ax.set_title(f'{h} — provisional reconstructed flame spectrum');ax.view_init(elev=27,azim=-135);fig.tight_layout();fig.savefig(figdir+f'/{h}_spectral_surface.png',dpi=100);plt.close(fig)
print('3D images complete',flush=True)
for f in sorted(glob.glob(R+'/curves_*.csv')):
 h=os.path.basename(f)[7:-4];d=pd.read_csv(f);t=d.time_s
 fig,axs=plt.subplots(3,2,figsize=(12,8),sharex=True)
 panels=[('oxygen_flow_spline','Oxygen-lance flow (m³/s)'),('oxygen_cumulative','Cumulative oxygen (m³)'),('CO2_flow','CO₂ instantaneous volume flow (m³/s)'),('CO_flow','CO instantaneous volume flow (m³/s)'),('CO2_cumulative','Cumulative CO₂ emission (m³)'),('CO_cumulative','Cumulative CO emission (m³)')]
 for ax,(col,title) in zip(axs.ravel(),panels):ax.plot(t,d[col],lw=.8);ax.set_title(title);ax.grid(alpha=.2)
 for ax in axs[-1,:]:ax.set_xlabel('Time (s)')
 fig.suptitle(h+' — process curves');fig.tight_layout();fig.savefig(figdir+f'/{h}_process_curves.png',dpi=100);plt.close(fig)
print('process images complete',flush=True)
