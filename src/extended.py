"""Extended biological-level inference and model sensitivity; generates numerical CSV only."""
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from scipy.optimize import least_squares
from .pipeline import DU_KEYS,CONC,DOSE,fit_hill,hill,fit_interaction,model_predict,calculate_bliss,summary

SEED=20261104

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan)
    good=np.flatnonzero(np.isfinite(p)); k=len(good)
    if not k:return out
    order=good[np.argsort(p[good])]; vals=p[order]*k/np.arange(1,k+1)
    out[order]=np.clip(np.minimum.accumulate(vals[::-1])[::-1],0,1)
    return out

def contrasts(du,hp,out):
    records=[]
    def add(meta, vals):
        x=np.asarray(vals,float)
        if len(x)!=3 or not np.isfinite(x).all():raise ValueError('Expected exactly three paired biological differences')
        sd=float(x.std(ddof=1))
        p=float(ttest_1samp(x,0).pvalue) if sd>1e-12 else (1. if abs(x.mean())<1e-12 else np.nan)
        signs=np.asarray(list(product((-1.,1.),repeat=3)))
        p_exact=float((np.abs(signs@x/3)>=abs(x.mean())-1e-10).mean())
        records.append({**meta,'n_bio':3,'difference_mean_pp':float(x.mean()),'difference_sd_pp':sd,'p_t_exploratory':p,'p_exact_signflip':p_exact,'all_same_direction':bool(np.all(x>0) or np.all(x<0)),**{f'exp{i+1}_pp':float(v) for i,v in enumerate(x)}})
    for key,g in du.groupby(DU_KEYS):
        p=g.pivot(index='bio_experiment',columns=[CONC,DOSE],values='response_pct')
        for dose in [0,10]:
            for c in [10,25,50,100]:
                add(dict(zip(DU_KEYS,key),concentration_uM=c,docetaxel_nM=dose,comparison='dose_vs_arm_baseline',family=f'{key}_DTX{dose}'),p[(c,dose)]-p[(0,dose)])
        for c in [10,25,50,100]:
            expected=p[(c,0)]*p[(0,10)]/p[(0,0)]
            add(dict(zip(DU_KEYS,key),concentration_uM=c,docetaxel_nM=10,comparison='observed_minus_bliss',family=f'{key}_bliss'),p[(c,10)]-expected)
            add(dict(zip(DU_KEYS,key),concentration_uM=c,docetaxel_nM=10,comparison='combo_minus_dtx',family=f'{key}_combo_dtx'),p[(c,10)]-p[(0,10)])
    h=hp.pivot(index='bio_experiment',columns='treatment_group',values='response_pct')
    for tr,ref in [('Apigenin','Control'),('Quercetin','Control'),('DTX+Apigenin','DTX'),('DTX+Quercetin','DTX')]:
        add(dict(compound=tr,time_h=24,assay='MTT',concentration_uM=50,docetaxel_nM=10 if tr.startswith('DTX+') else 0,comparison='hprec_treatment_vs_reference',family='HPrEC_four_prespecified'),h[tr]-h[ref])
    r=pd.DataFrame(records)
    r['q_BH_within_family']=np.nan;r['q_exact_BH_within_family']=np.nan
    for _,ix in r.groupby('family').groups.items():
        r.loc[ix,'q_BH_within_family']=bh(r.loc[ix,'p_t_exploratory'])
        r.loc[ix,'q_exact_BH_within_family']=bh(r.loc[ix,'p_exact_signflip'])
    r['q_BH_all_tests']=bh(r.p_t_exploratory)
    r.to_csv(out/'paired_biological_contrasts.csv',index=False)
    r.groupby(['comparison'],as_index=False).agg(n=('comparison','size'),min_exact_p=('p_exact_signflip','min'),same_direction_n=('all_same_direction','sum')).to_csv(out/'inference_sensitivity_summary.csv',index=False)

def lobo(du,out):
    rows=[];fitted=[]
    for key,g in du[du[DOSE].eq(0)].groupby(DU_KEYS):
        mat=g.pivot(index='bio_experiment',columns=CONC,values='response_pct').sort_index(axis=1)
        x=mat.columns.to_numpy(float)
        for held in mat.index:
            train=mat.drop(held).mean(axis=0).to_numpy(float)
            frame=pd.DataFrame({CONC:x,DOSE:0,'mean_pct':train})
            p=fit_hill(frame)
            if p is None:raise RuntimeError(f'Hill LOBO fit failed: {key}')
            pred=hill(x,p['bottom_pct'],p['ec50_uM'],p['hill_slope'],p['top_pct'])
            actual=mat.loc[held].to_numpy(float)
            rows.append(dict(zip(DU_KEYS,key),heldout_bio_experiment=int(held),mae_pp=float(np.mean(abs(pred-actual))),rmse_pp=float(np.sqrt(np.mean((pred-actual)**2))),train_ec50_uM=p['ec50_uM'],train_absolute_ic50_uM=p['absolute_ic50_uM']))
            for c,a,y in zip(x,actual,pred):fitted.append(dict(zip(DU_KEYS,key),heldout_bio_experiment=int(held),concentration_uM=float(c),observed_pct=float(a),predicted_pct=float(y)))
    pd.DataFrame(rows).to_csv(out/'lobo_hill_metrics.csv',index=False)
    pd.DataFrame(fitted).to_csv(out/'lobo_hill_predictions.csv',index=False)

def per_experiment_ic50(du,out):
    rows=[]
    for key,g in du[du[DOSE].eq(0)].groupby(DU_KEYS+['bio_experiment']):
        p=fit_hill(g.rename(columns={'response_pct':'mean_pct'}))
        if p is None:continue
        rows.append({**dict(zip(DU_KEYS+['bio_experiment'],key)),**p})
    pd.DataFrame(rows).to_csv(out/'hill_per_biological_experiment.csv',index=False)

def linear(x,y,t):
    x=np.asarray(x,float);y=np.asarray(y,float)
    order=np.argsort(x)
    return float(np.interp(t,x[order],y[order]))

def heldout_benchmarks(du,out):
    s=summary(du);rows=[]
    for key,g in s.groupby(DU_KEYS):
        for held in [10.,25.,50.]:
            train=g[g[CONC].ne(held)]; test=g[g[CONC].eq(held)]
            # Reproduce the original Phase 2 optimizer constraints and initial values.
            single_train=train[train[DOSE].eq(0)].sort_values(CONC)
            combo_train=train[train[DOSE].eq(10)].sort_values(CONC)
            ctrl0=float(single_train[single_train[CONC].eq(0)].mean_pct.iloc[0])
            doc0=float(combo_train[combo_train[CONC].eq(0)].mean_pct.iloc[0])
            xs=single_train[CONC].to_numpy(float); ys=single_train.mean_pct.to_numpy(float)
            def hh(x,z):return hill(x,z[0],z[1],z[2],ctrl0)
            z=least_squares(lambda v:hh(xs,v)-ys,[max(0,min(ys)*.75),35.,1.2],bounds=([0.,.1,.2],[min(150.,max(ctrl0,.1)),1000.,6.]),max_nfev=3000).x
            xx=combo_train[CONC].to_numpy(float)
            rr=combo_train.mean_pct.to_numpy(float)-hh(xx,z)*doc0/ctrl0
            mask=xx>0
            fit=least_squares(lambda v:np.r_[v[0]*xx[mask]/(v[1]+xx[mask])-rr[mask],.15*v[0]],
                [float(np.mean(rr[mask])),35.],bounds=([-100.,1.],[100.,500.]),max_nfev=2000)
            amp,k=fit.x
            single=train[train[DOSE].eq(0)]; combo=train[train[DOSE].eq(10)]
            actual=float(test[test[DOSE].eq(10)].mean_pct.iloc[0]); ctrl=float(single[single[CONC].eq(0)].mean_pct.iloc[0]);dtx=float(combo[combo[CONC].eq(0)].mean_pct.iloc[0])
            pb=float(hh(held,z)*doc0/ctrl0);pi=float(pb+amp*held/(k+held))
            preds={'hill_bliss':pb,'hill_interaction':pi,'direct_linear_combination':linear(combo[CONC],combo.mean_pct,held),'linear_bliss':linear(single[CONC],single.mean_pct,held)*dtx/ctrl}
            for model,pred in preds.items():rows.append({**dict(zip(DU_KEYS,key)),'heldout_uM':held,'model':model,'observed_pct':actual,'predicted_pct':pred,'error_pp':pred-actual})
    pred=pd.DataFrame(rows);pred.to_csv(out/'all_holdout_predictions.csv',index=False)
    metrics=pred.groupby('model',as_index=False).agg(n=('error_pp','size'),mae_pp=('error_pp',lambda x:float(np.mean(abs(x)))),rmse_pp=('error_pp',lambda x:float(np.sqrt(np.mean(x*x)))),bias_pp=('error_pp','mean'))
    metrics.to_csv(out/'all_holdout_metrics.csv',index=False)
    return metrics

def run(du_path,hp_path,out):
    out.mkdir(parents=True,exist_ok=True)
    du=pd.read_csv(du_path);hp=pd.read_csv(hp_path)
    contrasts(du,hp,out);lobo(du,out);per_experiment_ic50(du,out);biological_cluster_bootstrap(du,out)
    return heldout_benchmarks(du,out)

def biological_cluster_bootstrap(du,out,iterations=400):
    """Replicate original Phase 6 resampling order, fixed seed, full biological clusters."""
    from scipy.optimize import least_squares
    rng=np.random.default_rng(SEED);rows=[]
    for key,g in du[du[DOSE].eq(0)].groupby(DU_KEYS):
        mat=g.pivot(index='bio_experiment',columns=CONC,values='response_pct').sort_index(axis=1)
        doses=mat.columns.to_numpy(float)
        def fit(y):
            ctrl=float(y[doses==0][0]); upper=min(150.,max(ctrl,.1))
            res=least_squares(lambda z:hill(doses,z[0],z[1],z[2],ctrl)-y,
                [max(0,min(y)*.75),35,1.2],bounds=([0,.1,.2],[upper,1000,6]),max_nfev=3000)
            if not res.success:return None
            b,e,s=res.x; ic=np.nan
            if b<50<ctrl:
                ratio=(ctrl-b)/(50-b)-1
                if ratio>0:
                    v=e*ratio**(1/s)
                    if np.isfinite(v) and 0<=v<=100:ic=float(v)
            return float(e),ic
        samples=[]
        for _ in range(iterations):
            y=mat.iloc[rng.integers(0,len(mat),len(mat))].mean(axis=0).to_numpy(float)
            p=fit(y)
            if p is not None:samples.append(p)
        a=np.asarray(samples,float);defined=a[:,1][np.isfinite(a[:,1])]
        rows.append({**dict(zip(DU_KEYS,key)),'bootstrap_success':len(a),'ec50_q025':float(np.quantile(a[:,0],.025)),'ec50_q975':float(np.quantile(a[:,0],.975)),'ic50_defined_fraction':len(defined)/len(a),'ic50_q025_conditional':float(np.quantile(defined,.025)) if len(defined) else np.nan,'ic50_q975_conditional':float(np.quantile(defined,.975)) if len(defined) else np.nan})
    pd.DataFrame(rows).to_csv(out/'hill_biological_cluster_bootstrap.csv',index=False)
