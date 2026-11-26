from pathlib import Path
# Reuse identical model and plot functions, without rerunning the main sweep.
root=Path(__file__).resolve().parent
exec((root/'reproduce.py').read_text().split('for name,hs in ')[0],globals())
hs=[h for h in D if h!='6B03205']
results=[]
for scenario in ['official_0p5s','hypothetical_1s']:
 if scenario=='hypothetical_1s': D['6C02057']['time_s']*=2
 for lag in [0,48]:
  start=lag
  pp={g:fit(hs,lag,g,start) for g in G}
  if scenario=='hypothetical_1s':plot(hs,f'att2to10_H1_lag{lag:02d}_four_panel',lag,pp,start=start)
  for g in G:
   yall=[];pall=[];individual=[]
   for h in hs:
    t,x,y=get(h,lag,start);y=y[g];yall.extend(y);pall.extend(pred(g,x,pp[g]));individual.append(met(y,pred(g,x,fit([h],lag,g,start)))['r2'])
   results.append(dict(scenario=scenario,lag_s=lag,gas=g,pooled_r2=met(yall,pall)['r2'],per_heat_mean_r2=np.mean(individual)))
pd.DataFrame(results).to_csv(root/'time_scale_sensitivity.csv',index=False)
print(pd.DataFrame(results).to_string(index=False))
