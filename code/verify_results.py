"""Numerical checks separate from the main analysis workflow."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from reanalyse import fit
R=Path(__file__).resolve().parents[1]
d=pd.read_csv(R/'data/processed_panel.csv').dropna(subset=['df','dws','dd']).reset_index(drop=True)
f=fit(d,'df',['dws','dd']);d=f['d']
X=np.column_stack([np.ones(len(d)),d[['dws','dd']].to_numpy(),pd.get_dummies(d.iso_code,drop_first=True,dtype=float),pd.get_dummies(d.year,drop_first=True,dtype=float)])
y=d.df.to_numpy();b=np.linalg.lstsq(X,y,rcond=None)[0];e=y-X@b
assert np.allclose(b[1:3],f['b'][:,0],atol=1e-10)
gi=f['gi'];G=f['G'];bread=np.linalg.inv(X.T@X);score=np.zeros((G,X.shape[1]));np.add.at(score,gi,X*e[:,None]);v=f['correction']*bread@score.T@score@bread
assert np.allclose(v[1:3,1:3],f['covs'][0],atol=1e-10)
w=np.random.default_rng(73).choice([-1.,1.],size=G);ys=X@b+e*w[gi];bs=np.linalg.lstsq(X,ys,rcond=None)[0];es=ys-X@bs
a=f['xt']@f['bread'][:,0];shift=float(np.sum(a*e*w[gi]));assert abs((bs[1]-b[1])-shift)<1e-10
s=np.zeros(G);np.add.at(s,gi,a*es);direct=np.sqrt(f['correction']*sum(s*s))
Q=np.linalg.qr(X)[0];A=np.zeros((G,Q.shape[1]));U=np.zeros_like(A);sc=np.zeros(G)
np.add.at(A,gi,a[:,None]*Q);np.add.at(U,gi,e[:,None]*Q);np.add.at(sc,gi,a*e)
fastscore=sc*w-A@U.T@w;fast=np.sqrt(f['correction']*sum(fastscore*fastscore));assert abs(fast-direct)<1e-10
t=pd.read_csv(R/'results/T03_allocation.csv');aux=pd.read_csv(R/'results/S06_demand_auxiliary.csv').iloc[0].beta
fwl=[]
for outcome in ['df','dother','dm','dbalance']:
    sub=t[t.outcome.eq(outcome)];bc=sub[sub.model.eq('Conditional')&sub.term.eq('dws')].iloc[0].beta;delta=sub[sub.model.eq('Conditional')&sub.term.eq('dd')].iloc[0].beta;bu=sub[sub.model.eq('Unconditional')&sub.term.eq('dws')].iloc[0].beta
    err=abs(bu-bc-delta*aux);assert err<1e-10;fwl.append(err)
for term,target in [('dws',-1),('dd',1)]:assert abs(t[t.model.eq('Conditional')&t.term.eq(term)].beta.sum()-target)<1e-10
assert abs(t[t.model.eq('Unconditional')].beta.sum()+1)<1e-10
n=pd.read_csv(R/'results/T04_normalized.csv')
for group in n.group.unique():
    for term,target in [('zws',-1),('zd',1)]:assert abs(n[n.group.eq(group)&n.term.eq(term)].beta.sum()-target)<1e-9
assert len(d)==4726 and d.iso_code.nunique()==199
assert len(pd.read_csv(R/'results/S09_leave_one_out.csv'))==199
out=dict(status='PASS',direct_vs_absorbed_beta_max_error=float(abs(b[1:3]-f['b'][:,0]).max()),direct_vs_absorbed_cov_max_error=float(abs(v[1:3,1:3]-f['covs'][0]).max()),bootstrap_draw_coefficient_error=float(abs((bs[1]-b[1])-shift)),bootstrap_draw_se_error=float(abs(fast-direct)),fwl_max_error=float(max(fwl)),main_observations=len(d),main_systems=G,leave_one_out_models=199)
(R/'results/verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
