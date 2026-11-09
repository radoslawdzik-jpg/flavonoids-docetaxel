"""Reproduce downstream analysis from biological-experiment means; no figures."""
import argparse,json
from pathlib import Path
import pandas as pd
import numpy as np
from .pipeline import DU_KEYS,CONC,DOSE,summary,calculate_bliss,fit_hill,fit_interaction,holdout

def run(du,hprec,out):
    out.mkdir(parents=True,exist_ok=True)
    b=pd.read_csv(du)
    required=DU_KEYS+[CONC,DOSE,'bio_experiment','response_pct']
    if any(c not in b for c in required):raise ValueError('Missing required DU-145 biological-mean columns')
    keys=DU_KEYS+[CONC,DOSE,'bio_experiment']
    if b.duplicated(keys).any():raise ValueError('Duplicate experiment-condition rows')
    if set(b.bio_experiment.unique())!={1,2,3}:raise ValueError('Expected biological experiments 1, 2, 3')
    if b.response_pct.isna().any():raise ValueError('Missing responses')
    s=summary(b)
    b.to_csv(out/'du145_biological_means.csv',index=False)
    s.to_csv(out/'du145_condition_summary.csv',index=False)
    bliss=calculate_bliss(b)
    bliss.to_csv(out/'du145_bliss_biological.csv',index=False)
    if not bliss.empty:
        bliss.groupby(DU_KEYS+[CONC],as_index=False).agg(mean_deviation_pct=('bliss_deviation_pct','mean'),sd_deviation_pct=('bliss_deviation_pct','std'),n_bio=('bliss_deviation_pct','size')).to_csv(out/'du145_bliss_summary.csv',index=False)
    rows=[]
    for key,g in s.groupby(DU_KEYS):
        p=fit_hill(g)
        if p is None:continue
        ip=fit_interaction(g,p)
        row=dict(zip(DU_KEYS,key),**p)
        if ip:row.update(ip)
        rows.append(row)
    pd.DataFrame(rows).to_csv(out/'model_parameters.csv',index=False)
    pred=holdout(s)
    if not pred.empty:
        pred.to_csv(out/'concentration_holdout_predictions.csv',index=False)
        pred.groupby('model',as_index=False).agg(n=('error_pct','size'),mae_pct=('error_pct',lambda x:float(np.mean(np.abs(x)))),rmse_pct=('error_pct',lambda x:float(np.sqrt(np.mean(x*x)))),bias_pct=('error_pct','mean')).to_csv(out/'concentration_holdout_metrics.csv',index=False)
    if hprec:
        hp=pd.read_csv(hprec)
        req=['time_h','treatment_group','bio_experiment','response_pct']
        if any(c not in hp for c in req):raise ValueError('Missing HPrEC columns')
        if hp.duplicated(req[:-1]).any():raise ValueError('Duplicate HPrEC rows')
        hp.to_csv(out/'hprec_biological_means.csv',index=False)
        hp.groupby(['time_h','treatment_group'],as_index=False).agg(mean_pct=('response_pct','mean'),sd_pct=('response_pct','std'),n_bio=('response_pct','size')).to_csv(out/'hprec_condition_summary.csv',index=False)
    (out/'run_metadata.json').write_text(json.dumps({'input_level':'biological experiment means (four technical replicates averaged)','source_raw_data_included':False,'figure_generation':False,'biological_n':3,'note':'Downstream numerical rerun; not a verification of upstream raw absorbance normalization, full supplementary statistics or GEO analyses. Model outputs require manuscript concordance audit.'},indent=2)+'\n')
    print('Generated downstream outputs in',out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--du145',type=Path,default=Path('data/processed/du145_biological_means.csv'));p.add_argument('--hprec',type=Path,default=Path('data/processed/hprec_biological_means.csv'));p.add_argument('--out',type=Path,default=Path('outputs/reproduced'));a=p.parse_args();run(a.du145,a.hprec,a.out)
