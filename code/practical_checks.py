"""Additional practical-condition and exposure-support checks from the frozen panel.
Run after reanalyse.py; no download, source edits or imputation.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from reanalyse import fit
R=Path(__file__).resolve().parents[1]
d=pd.read_csv(R/'data/processed_panel.csv')
ys=['zf','zh','zn','zo','zm','zbalance']
c=d.dropna(subset=ys+['zws','zd']).copy()
rows=[]
# Baseline scale measured before the response period for all thresholds.
base=d[d.year.between(2000,2004)].groupby('iso_code').g.agg(['mean','count'])
for threshold in [0,.1,1.,5.]:
 ids=base.index[(base['mean']>=threshold)&(base['count']>=3)]
 f=fit(c[c.iso_code.isin(ids)],ys,['zws','zd'])
 rows.extend(dict(**x,condition=f'Baseline generation at least {threshold:g} TWh') for x in f['rows'])
for label,mask in [('Demand growth',c.dd>0),('Demand flat or declining',c.dd<=0)]:
 f=fit(c[mask],ys,['zws','zd'])
 rows.extend(dict(**x,condition=label) for x in f['rows'])
pd.DataFrame(rows).to_csv(R/'results/S18_practical_conditions.csv',index=False)
# Exposure variation for a partial slope: residualize exposure against FE and demand.
leverage=[]
for label,x,y,control in [('Absolute','dws','df','dd'),('Normalized','zws','zf','zd')]:
 f=fit(c,y,[x,control]);xx=f['xt'][:,0];z=f['xt'][:,1]
 v=xx-z*(z@xx)/(z@z)
 w=pd.DataFrame({'iso_code':f['d'].iso_code,'partial_ss':v*v})
 t=w.groupby('iso_code',as_index=False).partial_ss.sum();t['exposure_ss_share']=t.partial_ss/t.partial_ss.sum();t['specification']=label
 t=t.merge(d[['iso_code','country']].drop_duplicates(),on='iso_code')
 leverage.append(t.sort_values('exposure_ss_share',ascending=False))
pd.concat(leverage,ignore_index=True).to_csv(R/'results/S19_exposure_concentration.csv',index=False)
# Within-equation wind-versus-solar equality tests, using their covariance.
f=fit(d,['dcoal','dgas','doil','dh','dn','do','dm','dbalance'],['dw','ds','dd'])
tests=[]
for j,y in enumerate(f['ys']):
 v=np.array([1.,-1.,0.]);delta=float(v@f['b'][:,j]);se=float(np.sqrt(v@f['covs'][j]@v));crit=stats.t.ppf(.975,f['G']-1)
 tests.append(dict(outcome=y,wind_minus_solar=delta,se=se,lo=delta-crit*se,hi=delta+crit*se,p=2*stats.t.sf(abs(delta/se),f['G']-1),N=len(f['d']),systems=f['G']))
pd.DataFrame(tests).to_csv(R/'results/S20_technology_contrasts.csv',index=False)
# Recheck original-only results and exact allocation for each new condition.
t=pd.DataFrame(rows)
for name,sub in t.groupby('condition'):
 for term,target in [('zws',-1),('zd',1)]:
  assert abs(sub[sub.term.eq(term)].beta.sum()-target)<1e-9,name
out={'status':'PASS','practical_groups':t.condition.nunique(),'common_normalized_N':len(c),'common_normalized_systems':c.iso_code.nunique(),'baseline_scale_period':'2000–2004, at least three observed years','demand_regimes':'current demand change; descriptive realized-state comparisons, not exogenous regimes','technology_tests':'wind minus solar within one joint equation, full clustered covariance'}
(R/'results/practical_verification.json').write_text(json.dumps(out,indent=2))
print(t[t.term.eq('zws')&t.outcome.isin(['zf','zh','zm'])].to_string(index=False))
print(pd.concat(leverage).groupby('specification',sort=False).head(5).to_string(index=False))
print(pd.DataFrame(tests).to_string(index=False))
from reanalyse import wild
norm=leverage[1]
checks=[]
for label,ids in [('Exclude highest normalized leverage system',norm.head(1).iso_code),('Exclude five highest normalized leverage systems',norm.head(5).iso_code)]:
 f=fit(c[~c.iso_code.isin(ids)],ys,['zws','zd'])
 checks.extend(dict(**x,condition=label) for x in f['rows'])
 for j in [0,4]:
  z=wild(f,j=j)
  print('Additional wild interval',label,ys[j],z)
  pd.DataFrame([dict(outcome=ys[j],condition=label,**z)]).to_csv(R/'results'/('S22_wild_'+('top1' if len(ids)==1 else 'top5')+'_'+ys[j]+'.csv'),index=False)
 print('Normalized omitted systems',label,ids.tolist())
pd.DataFrame(checks).to_csv(R/'results/S21_normalized_influence.csv',index=False)
print(pd.DataFrame(checks).query("term=='zws' and outcome in ['zf','zm']").to_string(index=False))
# Source-derived wind/solar fuel component closure (never alter source values).
d['fossil_detail_residual']=d.f-d.coal-d.gas-d.oil
print('Fuel detail max discrepancy',d.fossil_detail_residual.abs().max())
# Independent direct-design check of representative additional comparisons.
errs=[]
for name,sub in [('Demand growth',c[c.dd>0]),('Demand flat or declining',c[c.dd<=0]),('Exclude highest normalized leverage system',c[c.iso_code.ne(norm.iloc[0].iso_code)])]:
 f=fit(sub,ys,['zws','zd']);s=f['d']
 X=np.column_stack([np.ones(len(s)),s[['zws','zd']],pd.get_dummies(s.iso_code,drop_first=True,dtype=float),pd.get_dummies(s.year,drop_first=True,dtype=float)])
 b=np.linalg.lstsq(X,s[ys].to_numpy(),rcond=None)[0][1:3]
 err=float(abs(b-f['b']).max());assert err<1e-9,(name,err);errs.append(err)
for name,sub in pd.DataFrame(checks).groupby('condition'):
 for term,target in [('zws',-1),('zd',1)]:assert abs(sub[sub.term.eq(term)].beta.sum()-target)<1e-9
out.update(direct_design_max_beta_error=max(errs),normalized_top_system=norm.iloc[0].country,normalized_top_exposure_share=float(norm.iloc[0].exposure_ss_share),fuel_detail_max_discrepancy=float(d.fossil_detail_residual.abs().max()))
(R/'results/practical_verification.json').write_text(json.dumps(out,indent=2))
