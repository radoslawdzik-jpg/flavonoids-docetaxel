#!/usr/bin/env python3
"""Data-only reproducibility pipeline for DU-145/HPrEC in-vitro responses.

No figures; no manuscript supplementary tables; no input observations distributed.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import ttest_rel

DU_KEYS=['compound','time_h','assay']
CONC='flavonoid_conc_uM'
RESP='viability_pct_vs_control'
DOSE='docetaxel_nM'
SEED=20261104

def check_columns(df, required, label):
    missing=sorted(set(required)-set(df.columns))
    if missing: raise ValueError(f'{label}: missing columns {missing}')

def load_du(path):
    d=pd.read_csv(path)
    check_columns(d,DU_KEYS+[CONC,DOSE,'replicate',RESP],'DU-145')
    d=d.copy()
    for col in [CONC,DOSE,'replicate',RESP,'time_h']:
        d[col]=pd.to_numeric(d[col],errors='raise')
    if d[[CONC,DOSE,'replicate',RESP]].isna().any().any():
        raise ValueError('Missing DU-145 dose, replicate or response')
    if (d[CONC]<0).any() or (d[DOSE]<0).any(): raise ValueError('Negative concentrations')
    if not d['replicate'].between(1,12).all() or not np.all(d['replicate']%1==0):
        raise ValueError('Expected DU-145 technical replicate IDs 1–12')
    d['bio_experiment']=((d['replicate'].astype(int)-1)//4)+1
    if set(d['bio_experiment'].unique())!={1,2,3}:
        raise ValueError('Expected three biological experiments')
    return d

def bio_means(d):
    cols=DU_KEYS+[CONC,DOSE,'bio_experiment']
    x=d.groupby(cols,as_index=False,dropna=False).agg(response_pct=(RESP,'mean'),n_technical=(RESP,'size'))
    if (x['n_technical']!=4).any():
        raise ValueError('Expected four technical replicates per biological experiment and condition; inspect duplicate or missing rows')
    return x

def summary(b):
    cols=DU_KEYS+[CONC,DOSE]
    s=b.groupby(cols,as_index=False).agg(mean_pct=('response_pct','mean'),sd_pct=('response_pct','std'),n_bio=('response_pct','size'))
    if (s.n_bio!=3).any(): raise ValueError('Expected n=3 biological experiments per condition')
    return s

def hill(x,bottom,ec50,slope,top):
    x=np.asarray(x,dtype=float)
    return bottom+(top-bottom)/(1+(x/np.maximum(ec50,1e-8))**slope)

def fit_hill(group):
    g=group[group[DOSE]==0].sort_values(CONC)
    g=g.groupby(CONC,as_index=False)['mean_pct'].mean()
    if len(g)<4 or not (g[CONC]==0).any(): return None
    top=float(g.loc[g[CONC]==0,'mean_pct'].mean())
    x=g[CONC].to_numpy(float); y=g.mean_pct.to_numpy(float)
    start=[max(0,min(float(np.min(y)),top)),50.,1.]
    def residual(p): return hill(x,*p,top)-y
    res=least_squares(residual,start,bounds=([0,1e-4,0.1],[max(150.,top+20),1e4,8]),max_nfev=10000)
    bottom,ec50,slope=map(float,res.x)
    # Absolute 50%-of-control threshold is only reported when observed within tested range.
    # This is separate from the Hill half-amplitude EC50.
    ic50=np.nan
    lo=float(np.min(x)); hi=float(np.max(x))
    if min(float(hill(lo,*res.x,top)),float(hill(hi,*res.x,top)))<=50<=max(float(hill(lo,*res.x,top)),float(hill(hi,*res.x,top))):
        from scipy.optimize import brentq
        ic50=float(brentq(lambda z:float(hill(z,*res.x,top))-50,lo,hi))
    return dict(top_pct=top,bottom_pct=bottom,ec50_uM=ec50,hill_slope=slope,absolute_ic50_uM=ic50,fit_rmse_pct=float(np.sqrt(np.mean(res.fun**2))))

def fit_interaction(g, pars, excluded=None):
    # Uses group-level means; optional exclusion of one flavonoid concentration.
    g=g.copy()
    if excluded is not None: g=g[g[CONC]!=excluded]
    if g.empty:return None
    ctrl=g[(g[CONC]==0)&(g[DOSE]==0)]
    dtx=g[(g[CONC]==0)&(g[DOSE]==10)]
    comb=g[(g[CONC]>0)&(g[DOSE]==10)]
    if ctrl.empty or dtx.empty or len(comb)<2:return None
    ctrl=float(ctrl.mean_pct.mean()); doc=float(dtx.mean_pct.mean())
    x=comb[CONC].to_numpy(float); obs=comb.mean_pct.to_numpy(float)
    baseline=hill(x,pars['bottom_pct'],pars['ec50_uM'],pars['hill_slope'],pars['top_pct'])*doc/ctrl
    def fun(p):
        A,K=p
        return np.r_[baseline+A*x/(K+x)-obs,0.15*A]
    res=least_squares(fun,[-10,35],bounds=([-150,0.01],[150,10000]),max_nfev=10000)
    return dict(amplitude_pct=float(res.x[0]),k_uM=float(res.x[1]),docetaxel_only_pct=doc,control_pct=ctrl,lambda_=0.15)

def model_predict(c,pars,ip):
    h=hill(c,pars['bottom_pct'],pars['ec50_uM'],pars['hill_slope'],pars['top_pct'])
    bliss=h*ip['docetaxel_only_pct']/ip['control_pct']
    interaction=bliss+ip['amplitude_pct']*c/(ip['k_uM']+c)
    return float(bliss),float(interaction)

def calculate_bliss(b):
    rows=[]
    for key,g in b.groupby(DU_KEYS+['bio_experiment']):
        control=g[(g[CONC]==0)&(g[DOSE]==0)]
        doc=g[(g[CONC]==0)&(g[DOSE]==10)]
        if len(control)!=1 or len(doc)!=1: raise ValueError(f'Missing or duplicated zero-dose control for {key}')
        ctrl=float(control.response_pct.iloc[0]); dtx=float(doc.response_pct.iloc[0])
        if ctrl<=0: raise ValueError('Control must be >0')
        for c in sorted(g[CONC].unique()):
            if c==0:continue
            single=g[(g[CONC]==c)&(g[DOSE]==0)]
            comb=g[(g[CONC]==c)&(g[DOSE]==10)]
            if len(single)!=1 or len(comb)!=1: continue
            obs=float(comb.response_pct.iloc[0]); exp=float(single.response_pct.iloc[0])*dtx/ctrl
            rows.append(dict(zip(DU_KEYS+['bio_experiment'],key),**{CONC:float(c)},observed_pct=obs,bliss_expected_pct=exp,bliss_deviation_pct=obs-exp))
    return pd.DataFrame(rows)

def holdout(s,holdouts=(10.,25.,50.)):
    rows=[]
    for key,g in s.groupby(DU_KEYS):
        for omitted in holdouts:
            if omitted not in set(g[CONC]):continue
            train=g[g[CONC]!=omitted]
            test=g[(g[CONC]==omitted)&(g[DOSE]==10)]
            if len(test)!=1:continue
            pars=fit_hill(train)
            if pars is None:continue
            ip=fit_interaction(train,pars)
            if ip is None:continue
            pred_b,pred_i=model_predict(omitted,pars,ip)
            obs=float(test.mean_pct.iloc[0])
            for model,pred in [('hill_bliss',pred_b),('hill_interaction',pred_i)]:
                rows.append(dict(zip(DU_KEYS,key),heldout_uM=omitted,model=model,observed_pct=obs,predicted_pct=pred,error_pct=pred-obs))
    return pd.DataFrame(rows)

def run(du_path,hprec_path,out):
    out.mkdir(parents=True,exist_ok=True)
    du=load_du(du_path); b=bio_means(du); s=summary(b)
    b.to_csv(out/'du145_biological_means.csv',index=False)
    s.to_csv(out/'du145_condition_summary.csv',index=False)
    bliss=calculate_bliss(b)
    bliss.to_csv(out/'du145_bliss_biological.csv',index=False)
    if not bliss.empty:
        bliss.groupby(DU_KEYS+[CONC],as_index=False).agg(mean_deviation_pct=('bliss_deviation_pct','mean'),sd_deviation_pct=('bliss_deviation_pct','std'),n_bio=('bliss_deviation_pct','size')).to_csv(out/'du145_bliss_summary.csv',index=False)
    models=[]
    for key,g in s.groupby(DU_KEYS):
        p=fit_hill(g)
        if p is None:continue
        ip=fit_interaction(g,p)
        row=dict(zip(DU_KEYS,key),**p)
        if ip:row.update(ip)
        models.append(row)
    pd.DataFrame(models).to_csv(out/'model_parameters.csv',index=False)
    pred=holdout(s)
    if not pred.empty:
        pred.to_csv(out/'concentration_holdout_predictions.csv',index=False)
        metrics=pred.groupby('model',as_index=False).agg(n=('error_pct','size'),mae_pct=('error_pct',lambda z:float(np.mean(np.abs(z)))),rmse_pct=('error_pct',lambda z:float(np.sqrt(np.mean(z**2)))),bias_pct=('error_pct','mean'))
        metrics.to_csv(out/'concentration_holdout_metrics.csv',index=False)
    if hprec_path:
        hp=pd.read_csv(hprec_path)
        check_columns(hp,['time_h','treatment_group','replicate',RESP],'HPrEC')
        hp['replicate']=pd.to_numeric(hp['replicate'],errors='raise')
        hp['bio_experiment']=((hp['replicate'].astype(int)-1)//4)+1
        if not hp['replicate'].between(1,12).all():raise ValueError('HPrEC replicate IDs must be 1–12')
        hp_b=hp.groupby(['time_h','treatment_group','bio_experiment'],as_index=False).agg(response_pct=(RESP,'mean'),n_technical=(RESP,'size'))
        if (hp_b.n_technical!=4).any():raise ValueError('HPrEC expected four technical replicates per biological experiment and group')
        hp_b.to_csv(out/'hprec_biological_means.csv',index=False)
        hp_b.groupby(['time_h','treatment_group'],as_index=False).agg(mean_pct=('response_pct','mean'),sd_pct=('response_pct','std'),n_bio=('response_pct','size')).to_csv(out/'hprec_condition_summary.csv',index=False)
    metadata={'analysis':'DU-145/HPrEC empirical treatment-response analysis','random_seed':SEED,'biological_n':3,'technical_per_bio':4,'input_data_included':False,'external_validation':False,'figure_generation':False,'bliss_synergy_claim':False,'note':'Internal concentration holdout is not external validation. Hill EC50 differs from absolute 50%-control IC50. Model specifications may differ from frozen manuscript publication fits.'}
    (out/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Analysis completed: {out}')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--du145',type=Path,required=True,help='Local authorised normalized DU-145 CSV')
    p.add_argument('--hprec',type=Path,help='Optional local authorised normalized HPrEC CSV')
    p.add_argument('--out',type=Path,default=Path('outputs'))
    a=p.parse_args();run(a.du145,a.hprec,a.out)
if __name__=='__main__':main()
