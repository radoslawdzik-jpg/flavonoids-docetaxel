#!/usr/bin/env python3
"""Run every available numerical workflow; figures are never created."""
import argparse, json, hashlib, sys, platform
from pathlib import Path
import numpy as np,pandas as pd,scipy
from src.run_processed import run as core
from src.extended import run as extended
from src.transcriptomics import main as geo_main

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--with-geo',action='store_true',help='Run GEO if external inputs are available; see README');p.add_argument('--out',type=Path,default=Path('outputs'));a=p.parse_args()
    du=Path('data/processed/du145_biological_means.csv');hp=Path('data/processed/hprec_biological_means.csv')
    if not du.exists() or not hp.exists():raise FileNotFoundError('Processed input files not found; run from repository root')
    out=a.out;out.mkdir(parents=True,exist_ok=True)
    core(du,hp,out/'cellular');metrics=extended(du,hp,out/'cellular')
    manifest={'input_sha256':{str(z):sha(z) for z in [du,hp]},'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'geo_executed':False,'figure_generation':False,'reproduction_level':'downstream numerical reanalysis, manuscript concordance not certified'}
    if a.with_geo:
        from src.transcriptomics import run
        ext=Path('data/external');args=[ext/'GSE83654_series_matrix.txt.gz',ext/'GPL13607-20416.txt',ext/'h.all.v2025.1.Hs.symbols.gmt',ext/'c2.cp.reactome.v2025.1.Hs.symbols.txt']
        if not all(z.exists() for z in args):raise FileNotFoundError('GEO requested, but one or more licensed external inputs are missing; see README')
        run(*args,out/'transcriptomics');manifest['geo_executed']=True
        manifest['external_sha256']={z.name:sha(z) for z in args}
    # Manuscript text benchmark (rounded to 2 decimals); not a full-table concordance test.
    expected={'hill_bliss':(5.54,6.69),'hill_interaction':(1.54,2.32),'direct_linear_combination':(2.26,2.89)}
    checks=[]
    for model,(mae,rmse) in expected.items():
        r=metrics.set_index('model').loc[model]
        checks.append({'metric':model,'manuscript_MAE_2dp':mae,'computed_MAE_2dp':round(float(r.mae_pp),2),'MAE_match':round(float(r.mae_pp),2)==mae,'manuscript_RMSE_2dp':rmse,'computed_RMSE_2dp':round(float(r.rmse_pp),2),'RMSE_match':round(float(r.rmse_pp),2)==rmse})
    pd.DataFrame(checks).to_csv(out/'holdout_manuscript_concordance.csv',index=False)
    (out/'reproduction_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Holdout metrics:');print(metrics.to_string(index=False))
    print('Reproduction status: numerical workflows ran; see docs/REPRODUCIBILITY_AUDIT.md for limitations.')
if __name__=='__main__':main()
