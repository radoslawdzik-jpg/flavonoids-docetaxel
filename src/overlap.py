"""Exploratory drug–gene overlap with GEO docetaxel time-course; no causal inference."""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

DRUGS=('docetaxel','apigenin','quercetin')

def bh(p):
    a=np.asarray(p,float); order=np.argsort(a); out=np.empty(len(a)); v=np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];out[order]=np.minimum(v,1);return out

def gmt(path):
    result={}
    for line in Path(path).read_text().splitlines():
        a=line.rstrip('\n').split('\t')
        if len(a)>=3:result[a[0]]=set(a[2:])
    return result

def run(associations,expression,hallmark,reactome,out,seed=42,permutations=1000):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    a=pd.read_csv(associations,dtype=str).fillna(''); required={'drug','gene_symbol','relation_type','source','evidence_id'}
    if not required.issubset(a.columns):raise ValueError('Association file missing: '+str(sorted(required-set(a.columns))))
    a['drug']=a.drug.str.strip().str.lower();a['gene_symbol']=a.gene_symbol.str.strip().str.upper()
    a=a[a.drug.isin(DRUGS)&a.gene_symbol.ne('')].drop_duplicates()
    if a.empty:raise ValueError('No real associations supplied. Do not use example/synthetic data for scientific results.')
    x=pd.read_csv(expression);x['gene_symbol']=x.gene_symbol.str.upper();x=x.drop_duplicates('gene_symbol').set_index('gene_symbol')
    times=[c for c in x if c.endswith('h_mean')]
    if not times:raise ValueError('No time columns (e.g. 24h_mean)')
    universe=set(x.index); sets={d:set(a.loc[a.drug.eq(d),'gene_symbol'])&universe for d in DRUGS}
    provenance=a.groupby(['drug','relation_type','source'],dropna=False).size().reset_index(name='records');provenance.to_csv(out/'association_provenance.csv',index=False)
    membership=pd.DataFrame({'gene_symbol':sorted(set.union(*(sets.values())))});
    for d in DRUGS:membership[d]=membership.gene_symbol.isin(sets[d]);
    membership.to_csv(out/'gene_overlap_membership.csv',index=False)
    from itertools import combinations
    overlaps=[]
    for k in (2,3):
        for combo in combinations(DRUGS,k):
            inter=set.intersection(*(sets[d] for d in combo)); n=len(universe); na=len(sets[combo[0]]); nb=len(sets[combo[1]])
            # Pairwise hypergeometric only; triple intersection descriptive.
            p=float(hypergeom.sf(len(inter)-1,n,na,nb)) if k==2 else np.nan
            overlaps.append({'comparison':' & '.join(combo),'n_overlap':len(inter),'genes':';'.join(sorted(inter)),'hypergeom_p':p,'background_n':n})
    ov=pd.DataFrame(overlaps); mask=ov.hypergeom_p.notna();ov.loc[mask,'q_BH']=bh(ov.loc[mask,'hypergeom_p']);ov.to_csv(out/'drug_gene_overlap.csv',index=False)
    pathways={**{'Hallmark:'+k:v for k,v in gmt(hallmark).items()},**{'Reactome:'+k:v for k,v in gmt(reactome).items()}}
    rows=[]
    for drug,genes in sets.items():
        for pathway,members in pathways.items():
            eligible=members&universe
            if not (15<=len(eligible)<=500) or not genes:continue
            intersection=genes&eligible
            p=hypergeom.sf(len(intersection)-1,len(universe),len(eligible),len(genes))
            rows.append({'drug':drug,'pathway':pathway,'drug_gene_n':len(genes),'pathway_n':len(eligible),'overlap_n':len(intersection),'overlap_genes':';'.join(sorted(intersection)),'p':float(p)})
    e=pd.DataFrame(rows)
    if not e.empty:
        # BH correction per drug across all tested gene sets
        e['q_BH']=np.nan
        for drug,idx in e.groupby('drug').groups.items():e.loc[idx,'q_BH']=bh(e.loc[idx,'p'])
    e.to_csv(out/'drug_pathway_enrichment.csv',index=False)
    # Descriptive overlap with time-course expression, without DE claims.
    z=[]
    for drug,genes in sets.items():
        for t in times:
            s=x.loc[sorted(genes),t].dropna().astype(float)
            if len(s):z.append({'drug':drug,'time':t,'n_genes':len(s),'median_signal':float(s.median()),'fraction_positive':float((s>0).mean()),'note':'Descriptive GEO signal; not a validated differential expression estimate'})
    pd.DataFrame(z).to_csv(out/'drug_gene_timecourse_descriptive.csv',index=False)
    # Pair drug-associated pathways with observed GEO pathway rank effects, when available.
    pd.DataFrame({'metric':['association_rows','background_genes','permutations_requested','seed'],'value':[len(a),len(universe),permutations,seed]}).to_csv(out/'run_summary.csv',index=False)
    (out/'interpretation_limits.txt').write_text('GSE83654 contains docetaxel-only two-colour technical dye swaps; it cannot establish flavonoid or combination-induced gene regulation. Associations differ in evidence type. Gene-set tests are exploratory and gene/pathway dependence is not resolved. No causal mechanism, synergy, or validated differential expression is inferred.\n')
    return {'association_rows':len(a),'background_genes':len(universe),'drug_counts':{d:len(s) for d,s in sets.items()},'pathway_tests':len(e)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--associations',required=True);p.add_argument('--expression',required=True);p.add_argument('--hallmark',required=True);p.add_argument('--reactome',required=True);p.add_argument('--out',default='outputs/overlap');args=p.parse_args();print(json.dumps(run(args.associations,args.expression,args.hallmark,args.reactome,args.out),indent=2))
if __name__=='__main__':main()
