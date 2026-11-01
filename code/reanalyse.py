"""Energy Reports EGYR-D-26-03711: revised descriptive accounting analysis.

Run: OPENBLAS_NUM_THREADS=1 python code/reanalyse.py --source data/owid-energy-data.xlsx
Requires numpy, pandas, scipy, openpyxl, matplotlib. No source values are imputed.
All regressions use explicit system/year fixed effects, CR1 cluster covariance,
and t(G-1) reference inference. Wild intervals use 9,999 unrestricted Rademacher
cluster bootstrap-t draws, with fixed design and re-estimated fixed effects.
"""
from pathlib import Path
import argparse, hashlib, json, shutil
import numpy as np
import pandas as pd
from scipy import stats, linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
for p in ['results','figures','data']: (ROOT/p).mkdir(exist_ok=True)
TABLES={}; META={}; FITS={}

def fit(data, ys, xs, trends=False):
    if isinstance(ys,str): ys=[ys]
    d=data.dropna(subset=ys+xs).copy().reset_index(drop=True)
    z=pd.concat([pd.get_dummies(d.iso_code,dtype=float),pd.get_dummies(d.year,drop_first=True,dtype=float)],axis=1).to_numpy()
    if trends:
        z=np.column_stack([z,pd.get_dummies(d.iso_code,dtype=float).to_numpy()*(d.year.to_numpy()-2012)[:,None]])
    q,rr,piv=linalg.qr(z,mode='economic',pivoting=True)
    rank=int(np.sum(abs(np.diag(rr))>abs(rr[0,0])*1e-10));q=q[:,:rank]
    x=d[xs].to_numpy(float); y=d[ys].to_numpy(float)
    xt=x-q@(q.T@x);yt=y-q@(q.T@y)
    bread=np.linalg.inv(xt.T@xt);b=bread@xt.T@yt;u=yt-xt@b
    groups,gi=np.unique(d.iso_code,return_inverse=True);G=len(groups);n=len(d)
    correction=G/(G-1)*(n-1)/(n-rank-len(xs))
    rows=[];covs=[]
    for j,yname in enumerate(ys):
        s=np.zeros((G,len(xs)));np.add.at(s,gi,xt*u[:,j,None])
        v=correction*bread@s.T@s@bread;covs.append(v)
        se=np.sqrt(np.maximum(np.diag(v),0));t=b[:,j]/np.maximum(se,1e-300);crit=stats.t.ppf(.975,G-1)
        for k,xname in enumerate(xs):
            rows.append(dict(outcome=yname,term=xname,beta=b[k,j],se=se[k],lo=b[k,j]-crit*se[k],hi=b[k,j]+crit*se[k],p=2*stats.t.sf(abs(t[k]),G-1),N=n,systems=G))
    return dict(rows=rows,d=d,q=q,xt=xt,bread=bread,b=b,u=u,gi=gi,G=G,correction=correction,covs=covs,xs=xs,ys=ys)

def add(name, f, **tags):
    rows=[dict(**r,**tags) for r in f['rows']];TABLES.setdefault(name,[]).extend(rows);return f

def wild(f,j=0,k=0,B=9999,seed=20261101):
    # a is the full-design influence vector for the selected slope. Projection
    # via an orthonormal basis re-estimates all FE in every bootstrap draw.
    a=f['xt']@f['bread'][:,k];u=f['u'][:,j];gi=f['gi'];G=f['G']
    qx=linalg.orth(f['xt']);Q=np.column_stack([f['q'],qx])
    A=np.zeros((G,Q.shape[1]));U=np.zeros((G,Q.shape[1]));s=np.zeros(G)
    np.add.at(A,gi,a[:,None]*Q);np.add.at(U,gi,u[:,None]*Q);np.add.at(s,gi,a*u)
    T=A@U.T; rng=np.random.default_rng(seed);ts=[]
    for start in range(0,B,1000):
        w=rng.choice([-1.,1.],size=(G,min(1000,B-start)))
        delta=s@w;sstar=s[:,None]*w-T@w
        sest=np.sqrt(f['correction']*np.sum(sstar*sstar,axis=0))
        ts.extend((delta/np.maximum(sest,1e-300)).tolist())
    ts=np.asarray(ts);b=f['b'][k,j];se=np.sqrt(f['covs'][j][k,k]);lo,hi=np.quantile(ts,[.025,.975])
    return dict(beta=b,lo=b-hi*se,hi=b-lo*se,p=(1+np.sum(abs(ts)>=abs(b/se)))/(B+1),draws=B,seed=seed,N=len(f['d']),systems=G)

def export(name,df):
    if not isinstance(df,pd.DataFrame):df=pd.DataFrame(df)
    df.to_csv(ROOT/'results'/f'{name}.csv',index=False)
    return df

def main(source):
    raw=pd.read_excel(source,sheet_name='data');metadata=pd.read_excel(source,sheet_name='metadata')
    META['source_sha256']=hashlib.sha256(Path(source).read_bytes()).hexdigest()
    r=raw[raw.iso_code.fillna('').str.len().eq(3)&raw.year.between(2000,2024)].copy()
    core=['wind_electricity','solar_electricity','fossil_electricity','electricity_demand']
    coverage=r[core].notna().all(axis=1).groupby(r.iso_code).sum()
    oldcoverage=r[core+['carbon_intensity_elec']].notna().all(axis=1).groupby(r.iso_code).sum()
    candidates=coverage[coverage>=20].index;oldelig=oldcoverage[oldcoverage>=20].index
    positive=r.groupby('iso_code').electricity_generation.max()
    elig=candidates.intersection(positive[positive>0].index)
    META['added_after_dropping_carbon_filter']=list(candidates.difference(oldelig))
    META['excluded_no_positive_generation']=list(candidates.difference(elig))
    r=r[r.iso_code.isin(elig)].sort_values(['iso_code','year']).copy()
    short={'f':'fossil_electricity','w':'wind_electricity','s':'solar_electricity','h':'hydro_electricity','n':'nuclear_electricity','o':'other_renewable_electricity','m':'net_elec_imports','d':'electricity_demand','g':'electricity_generation','coal':'coal_electricity','gas':'gas_electricity','oil':'oil_electricity'}
    for k,c in short.items():r[k]=r[c]
    r['ws']=r.w+r.s;r['other']=r.h+r.n+r.o;r['balance']=r.d-r.f-r.ws-r.other-r.m
    r['generation_residual']=r.g-r.f-r.ws-r.other;r['demand_residual']=r.d-r.g-r.m
    means=r.groupby('iso_code').g.mean()
    # Rank groups, with surplus systems assigned to earlier groups: 50/50/50/49.
    ordered=means.rename('mean').reset_index().sort_values(['mean','iso_code']).iso_code.to_numpy()
    quart={code:f'Q{j+1}' for j,ids in enumerate(np.array_split(ordered,4)) for code in ids}
    r['quartile']=r.iso_code.map(quart)
    for c in [*short,'ws','other','balance']:
        r['d'+c]=r.groupby('iso_code')[c].diff().where(r.groupby('iso_code').year.diff().eq(1))
    r['lagg']=r.groupby('iso_code').g.shift();r.loc[r.lagg<=0,'lagg']=np.nan
    for c in [*short,'ws','other','balance']:r['z'+c]=100*r['d'+c]/r.lagg
    r['minus_dd']=-r.dd;r['minus_zd']=-r.zd
    for lag in [1,2]:r['lag'+str(lag)]=r.groupby('iso_code').dws.shift(lag).where(r.groupby('iso_code').year.diff(lag).eq(lag))
    r['lead1']=r.groupby('iso_code').dws.shift(-1).where(r.groupby('iso_code').year.diff(-1).eq(-1))
    panel=r.dropna(subset=['df','dws','dd']).copy()
    common=panel.dropna(subset=['dother','dm','dbalance']).copy()
    balanced=r[core+['carbon_intensity_elec']].notna().all(axis=1).groupby(r.iso_code).sum();balids=balanced[balanced==25].index
    META.update(eligible=len(elig),oldeligible=len(oldelig),mainN=len(panel),commonN=len(common),commonG=common.iso_code.nunique(),balancedG=len(balids),kosovo_codes=raw.loc[raw.country.eq('Kosovo'),'iso_code'].dropna().unique().tolist())
    export('S01_sample_flow',[dict(stage='Three-character source identifiers, 2000–2024',systems=coverage.size,rows=len(raw[raw.iso_code.fillna('').str.len().eq(3)&raw.year.between(2000,2024)])),dict(stage='At least 20 joint core years, no carbon filter',systems=len(candidates),rows=int(raw.iso_code.isin(candidates).mul(raw.year.between(2000,2024)).sum())),dict(stage='Exclude systems with no positive recorded generation',systems=len(elig),rows=len(r)),dict(stage='Old rule including carbon intensity for comparison',systems=len(oldelig),rows=int(r.iso_code.isin(oldelig).sum())),dict(stage='Consecutive-year main differences',systems=panel.iso_code.nunique(),rows=len(panel)),dict(stage='Common accounting sample',systems=common.iso_code.nunique(),rows=len(common)),dict(stage='Common normalized accounting sample',systems=common.dropna(subset=['lagg']).iso_code.nunique(),rows=len(common.dropna(subset=['lagg']))),dict(stage='Balanced descriptive panel including carbon intensity',systems=len(balids),rows=len(balids)*25)])
    export('S02_system_coverage',r.groupby(['country','iso_code','quartile']).agg(years=('year','size'),mean_generation=('g','mean'),first_year=('year','min'),last_year=('year','max')).reset_index())
    metadata[metadata.column.isin(list(short.values())+['carbon_intensity_elec'])].to_csv(ROOT/'results'/'S03_source_metadata.csv',index=False)
    for c in ['generation_residual','demand_residual','balance']:
        x=r[c].dropna();META[c]=dict(n=len(x),max_abs=float(abs(x).max()),median_abs=float(abs(x).median()))
    export('S04_balance_residuals',r[['country','iso_code','year','generation_residual','demand_residual','balance']])
    add('T01_main',fit(panel,'df',['dws','dd']),model='Conditional')
    add('T01_main',fit(panel,'df',['dws']),model='Unconditional')
    for name,sub in [('Overall',panel),('Exclude China',panel[panel.iso_code.ne('CHN')]),('Exclude five largest',panel[~panel.iso_code.isin(means.nlargest(5).index)])]+[(q,panel[panel.quartile.eq(q)]) for q in ['Q1','Q2','Q3','Q4']]:
        f=add('T02_scale',fit(sub,'df',['dws','dd']),group=name)
        if name in ['Overall','Exclude China','Exclude five largest']:
            print('Wild bootstrap',name,flush=True);TABLES.setdefault('S05_wild_bootstrap',[]).append(dict(group=name,**wild(f)))
    META['largest_five']=r[r.iso_code.isin(means.nlargest(5).index)][['country','iso_code']].drop_duplicates().to_dict('records')
    f=fit(panel,'zf',['zws','zd']);add('T02_scale',f,group='Normalized overall')
    TABLES['S05_wild_bootstrap'].append(dict(group='Normalized overall',**wild(f)))
    for cond in [True,False]:
        ys=['df','dother','dm','dbalance']+([] if cond else ['minus_dd'])
        add('T03_allocation',fit(common,ys,['dws']+(['dd'] if cond else [])),model='Conditional' if cond else 'Unconditional')
    add('S06_demand_auxiliary',fit(common,'dd',['dws']))
    for name,sub in [('Overall',common)]+[(q,common[common.quartile.eq(q)]) for q in ['Q1','Q2','Q3','Q4']]:
        add('T04_normalized',fit(sub,['zf','zh','zn','zo','zm','zbalance'],['zws','zd']),group=name)
        add('S07_normalized_unconditional',fit(sub,['zf','zh','zn','zo','zm','zbalance','minus_zd'],['zws']),group=name)
    pair=panel.dropna(subset=['lag1']).copy()
    for c in ['dws','lag1']:pair[c+'_within']=pair[c]-pair.groupby('iso_code')[c].transform('mean')
    META['pooled_within_autocorrelation']=pair.dws_within.corr(pair.lag1_within)
    cors=pair.groupby('iso_code').apply(lambda d:d.dws.corr(d.lag1) if d.dws.std()>1e-12 and d.lag1.std()>1e-12 else np.nan,include_groups=False)
    META['median_system_autocorrelation']=float(cors.median());META['autocorr_systems']=int(cors.notna().sum())
    export('S08_autocorrelation',cors.rename('lag1_correlation').reset_index())
    for cond in [True,False]:
        for name,xs in [('Current',['dws']),('Lag only',['lag1']),('Lead only',['lead1']),('Lead current lag jointly',['lead1','dws','lag1']),('Distributed lag 0 to 2',['dws','lag1','lag2'])]:
            f=add('T05_timing',fit(panel,'df',xs+(['dd'] if cond else [])),model=name,demand_control=cond)
            if name=='Distributed lag 0 to 2':
                c=np.array([1.,1.,1.]+([0.] if cond else []));b=float(c@f['b'][:,0]);se=float(np.sqrt(c@f['covs'][0]@c));crit=stats.t.ppf(.975,f['G']-1)
                TABLES['T05_timing'].append(dict(outcome='df',term='cumulative_0_2',beta=b,se=se,lo=b-crit*se,hi=b+crit*se,p=2*stats.t.sf(abs(b/se),f['G']-1),N=len(f['d']),systems=f['G'],model=name,demand_control=cond))
    for code in sorted(panel.iso_code.unique()):
        f=fit(panel[panel.iso_code.ne(code)],'df',['dws','dd']);rr=f['rows'][0];TABLES.setdefault('S09_leave_one_out',[]).append(dict(excluded_code=code,excluded_system=panel.loc[panel.iso_code.eq(code),'country'].iloc[0],**rr))
    for q,sub in panel.groupby('quartile'):
        TABLES.setdefault('S10_zero_exposure',[]).append(dict(quartile=q,N=len(sub),systems=sub.iso_code.nunique(),zero_share=float(sub.dws.eq(0).mean()),abs_change_below_001_share=float(sub.dws.abs().lt(.01).mean()),median_abs_change=float(sub.dws.abs().median())))
    maxshare=(100*r.ws/r.g.where(r.g>0)).groupby(r.iso_code).max()
    for threshold in [.1,1.]:add('S11_sensitivity',fit(panel[panel.iso_code.isin(maxshare[maxshare>=threshold].index)],'df',['dws','dd']),model=f'Max observed WS share at least {threshold}%')
    add('S11_sensitivity',fit(panel[panel.year.lt(2024)],'df',['dws','dd']),model='Exclude 2024')
    add('S11_sensitivity',fit(panel,'df',['dws','dd'],trends=True),model='System trends in differences')
    add('T06_technology',fit(panel,['dcoal','dgas','doil','dh','dn','do','dm','dbalance'],['dw','ds','dd']),model='Wind and solar jointly')
    base=r[r.year.between(2000,2004)].copy();base['fshare']=100*base.f/base.g.where(base.g>0);base['hshare']=100*base.h/base.g.where(base.g>0);base['importshare']=100*base.m.abs()/base.d.where(base.d>0)
    for var in ['fshare','hshare','importshare']:
        bs=base.groupby('iso_code')[var].agg(['mean','count']);bs=bs[bs['count']>=3];cut=bs['mean'].median();META[var+'_cut']=float(cut)
        for label,ids in [('Lower',bs.index[bs['mean']<cut]),('Higher',bs.index[bs['mean']>=cut])]:
            add('S12_composition_groups',fit(panel[panel.iso_code.isin(ids)],['zf','zh','zn','zo','zm','zbalance'],['zws','zd']),group=label,baseline=var,cut=cut)
    desc=[]
    for c in ['dws','df','dd','zws','zf']:
        x=panel[c].dropna();desc.append(dict(variable=c,N=len(x),mean=x.mean(),sd=x.std(),median=x.median(),q1=x.quantile(.25),q3=x.quantile(.75)))
    export('S13_descriptive',desc)
    br=r[r.iso_code.isin(balids)].copy();tra=[]
    for year,sub in br.groupby('year'):
        tra.append(dict(year=year,N=len(sub),median_ws=sub.wind_share_elec.add(sub.solar_share_elec).median(),q1_ws=sub.wind_share_elec.add(sub.solar_share_elec).quantile(.25),q3_ws=sub.wind_share_elec.add(sub.solar_share_elec).quantile(.75),median_f=sub.fossil_share_elec.median(),q1_f=sub.fossil_share_elec.quantile(.25),q3_f=sub.fossil_share_elec.quantile(.75),weighted_ws=100*sub.ws.sum()/sub.g.sum(),weighted_f=100*sub.f.sum()/sub.g.sum()))
    export('S14_trajectories',tra)
    endpoints=br[br.year.isin([2000,2024])].pivot(index=['country','iso_code'],columns='year',values=['fossil_share_elec','wind_share_elec','solar_share_elec','g','f','h','n','ws'])
    endpoints.columns=[f'{c}_{y}' for c,y in endpoints.columns];endpoints=endpoints.reset_index();endpoints['fossil_share_change']=endpoints.fossil_share_elec_2024-endpoints.fossil_share_elec_2000;endpoints['ws_share_change']=endpoints.wind_share_elec_2024+endpoints.solar_share_elec_2024-endpoints.wind_share_elec_2000-endpoints.solar_share_elec_2000
    export('S15_endpoint_changes',endpoints)
    extremes=endpoints[endpoints.fossil_share_change.abs()>65];export('S16_extreme_endpoint_systems',extremes)
    export('S17_extreme_system_annual_series',r[r.iso_code.isin(extremes.iso_code)][['country','iso_code','year','g','f','h','n','ws','fossil_share_elec','generation_residual','demand_residual']])
    add('S11_sensitivity',fit(panel[~panel.iso_code.isin(extremes.iso_code)],'df',['dws','dd']),model='Exclude extreme endpoint systems')
    r.to_csv(ROOT/'data'/'processed_panel.csv',index=False)
    for name,rows in TABLES.items():export(name,rows)
    (ROOT/'results'/'metadata.json').write_text(json.dumps(META,indent=2,ensure_ascii=False),encoding='utf8')
    make_figures()
    print(json.dumps(META,indent=2),flush=True)

def make_figures():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
    def save(fig,name):
        fig.savefig(ROOT/'figures'/(name+'.png'),dpi=300);fig.savefig(ROOT/'figures'/(name+'.pdf'));plt.close(fig)
    t=pd.read_csv(ROOT/'results'/'T02_scale.csv');t=t[t.term.eq('dws')]
    fig,ax=plt.subplots(figsize=(7.2,4.3));y=np.arange(len(t));ax.errorbar(t.beta,y,xerr=[t.beta-t.lo,t.hi-t.beta],fmt='o',color='#175b69',capsize=3);ax.set_yticks(y,t.group);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.8);ax.axvline(-1,color='gray',lw=.8,ls=':');ax.set_xlabel('Conditional fossil coefficient (TWh/TWh)');save(fig,'Figure_1_Scale')
    t=pd.read_csv(ROOT/'results'/'S14_trajectories.csv');fig,axes=plt.subplots(1,2,figsize=(8,3.3))
    for ax,c,l in zip(axes,['ws','f'],['Wind and solar share (%)','Fossil share (%)']):
        ax.fill_between(t.year,t['q1_'+c],t['q3_'+c],color='#d8e6e8');ax.plot(t.year,t['median_'+c],label='Unweighted median',color='#175b69');ax.plot(t.year,t['weighted_'+c],label='Generation-weighted',color='#9a4c20',ls='--');ax.set_ylabel(l);ax.set_xlabel('Year')
    axes[0].legend(fontsize=8);save(fig,'Figure_4_Trajectories')
    t=pd.read_csv(ROOT/'results'/'T04_normalized.csv');t=t[t.term.eq('zws')&t.outcome.ne('zbalance')];fig,ax=plt.subplots(figsize=(8,4));labels={'zf':'Fossil','zh':'Hydro','zn':'Nuclear','zo':'Other renewables','zm':'Net imports'}
    for k,(out,label) in enumerate(labels.items()):
        sub=t[t.outcome.eq(out)].set_index('group').loc[['Overall','Q1','Q2','Q3','Q4']];x=np.arange(5)+(k-2)*.13;ax.errorbar(x,sub.beta,yerr=[sub.beta-sub.lo,sub.hi-sub.beta],fmt='o',ms=4,capsize=2,label=label)
    ax.set_xticks(np.arange(5),['Overall','Q1','Q2','Q3','Q4']);ax.axhline(0,color='gray',lw=.7);ax.set_ylabel('Normalized conditional allocation coefficient');ax.legend(ncol=3,fontsize=8);save(fig,'Figure_3_Normalized_allocation')
    t=pd.read_csv(ROOT/'results'/'S09_leave_one_out.csv').sort_values('beta').reset_index(drop=True);fig,ax=plt.subplots(figsize=(8,4));ax.plot(np.arange(len(t)),t.beta,'o',ms=2,color='#175b69');ax.axhline(-.9539727163,color='gray',ls='--');
    for ix in list(t.index[:2])+list(t.index[-3:]):ax.annotate(t.loc[ix,'excluded_system'],(ix,t.loc[ix,'beta']),xytext=(0,12 if ix%2 else -16),textcoords='offset points',fontsize=8,ha='center')
    ax.set_xlabel(f'All {len(t)} leave-one-system-out estimates, sorted');ax.set_ylabel('Conditional fossil coefficient');save(fig,'Figure_2_Leave_one_out')
    t=pd.read_csv(ROOT/'results'/'S15_endpoint_changes.csv');fig,ax=plt.subplots(figsize=(7,4));ax.scatter(t.ws_share_change,t.fossil_share_change,s=18,color='#175b69',alpha=.6)
    for row in t[t.fossil_share_change.abs()>65].itertuples():ax.annotate(row.country,(row.ws_share_change,row.fossil_share_change),xytext=(8,4),textcoords='offset points',fontsize=8)
    ax.axhline(0,color='gray',lw=.7);ax.axvline(0,color='gray',lw=.7);ax.set_xlabel('Wind and solar share change, 2000–2024 (pp)');ax.set_ylabel('Fossil share change, 2000–2024 (pp)');save(fig,'Figure_S1_Endpoint_changes')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'data'/'owid-energy-data.xlsx');args=parser.parse_args();main(args.source)
