"""GSE83654 exploratory gene-set rank analysis. No source data are bundled."""
from pathlib import Path
import argparse, gzip, io, re, urllib.request
import numpy as np,pandas as pd
from scipy.stats import rankdata
TIMES={8:['GSM2211839','GSM2211840'],24:['GSM2211841','GSM2211842'],48:['GSM2211843','GSM2211844'],72:['GSM2211845','GSM2211846']}
URL_MATRIX='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83654/matrix/GSE83654_series_matrix.txt.gz'
URL_GPL='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13607/annot/GPL13607.annot.gz'

def download(url,path):
    if not path.exists():
        print('Downloading external data from:',url,flush=True)
        urllib.request.urlretrieve(url,path)

def read_matrix(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
        else:raise ValueError('GEO matrix table not found')
        header=next(f).strip().replace('"','').split('\t');rows=[]
        for line in f:
            if line.startswith('!series_matrix_table_end'):break
            rows.append(line.rstrip('\n').split('\t'))
    d=pd.DataFrame(rows,columns=header)
    for c in d.columns[1:]:d[c]=pd.to_numeric(d[c],errors='coerce')
    return d

def read_annotation(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',errors='replace') as f:
        lines=f.readlines()
    a=next((i for i,x in enumerate(lines) if x.startswith('!platform_table_begin')),None)
    if a is not None:
        b=next((i for i in range(a+1,len(lines)) if lines[i].startswith('!platform_table_end')),len(lines))
        data=''.join(lines[a+1:b])
    else:data=''.join(x for x in lines if not x.startswith(('!','#','^')))
    return pd.read_csv(io.StringIO(data),sep='\t',dtype=str,on_bad_lines='skip',low_memory=False)

def run(matrix,gpl,hallmark,reactome,out):
    out.mkdir(parents=True,exist_ok=True)
    m=read_matrix(matrix);a=read_annotation(gpl)
    if not {'ID','GeneName','ControlType'}.issubset(a.columns):
        raise ValueError('Expected GPL13607-20416 full platform table with ID, GeneName, ControlType. Supply original GPL13607-20416.txt; generic .annot.gz is not equivalent.')
    a=a[a.ControlType.eq('0')][['ID','GeneName']].rename(columns={'ID':'ID_REF','GeneName':'gene_symbol'})
    a.gene_symbol=a.gene_symbol.str.strip()
    a=a[a.gene_symbol.str.fullmatch(r'[A-Za-z][A-Za-z0-9-]{1,29}',na=False)]
    a=a[~a.gene_symbol.str.match(r'^(?:ENST\d+|A_\d+_P\d+|DCP_\d+|LINC?RNA:)',case=False,na=False)].drop_duplicates('ID_REF')
    probe=m[['ID_REF']].copy()
    for t,ids in TIMES.items():
        if not set(ids).issubset(m.columns):raise ValueError(f'Missing GSM columns for {t}h')
        probe[f'{t}h_mean']=m[ids].mean(axis=1)
        probe[f'{t}h_pair_absdiff']=(m[ids[0]]-m[ids[1]]).abs()
    merged=probe.merge(a,on='ID_REF',validate='one_to_one')
    gene=merged.groupby('gene_symbol')[[f'{t}h_mean' for t in TIMES]].median().dropna()
    gene.to_csv(out/'gene_timecourse_derived.csv',index_label='gene_symbol')
    sets=[]
    for collection,path in [('Hallmark',hallmark),('Reactome',reactome)]:
        for line in path.read_text().splitlines():
            z=line.split('\t')
            if len(z)>3:sets.append((collection,z[0],set(z[2:])))
    genes=gene.index.tolist(); idx={v:i for i,v in enumerate(genes)};n=len(genes)
    eligible=[]
    for coll,name,s in sets:
        ids=np.asarray([idx[v] for v in s if v in idx],int)
        if 15<=len(ids)<=500:eligible.append((coll,name,ids,len(s)))
    rows=[]
    for t in TIMES:
        vals=gene[f'{t}h_mean'].to_numpy(float);ranks=rankdata(vals,method='average')
        for coll,name,ids,size in eligible:
            k=len(ids);auc=(ranks[ids].sum()-k*(k+1)/2)/(k*(n-k))
            rows.append(dict(time_h=t,collection=coll,pathway=name,overlap_n=k,set_n=size,rank_auc=float(auc),rank_biserial=float(2*auc-1),median_response=float(np.median(vals[ids])),mean_response=float(np.mean(vals[ids])),fraction_positive=float(np.mean(vals[ids]>0))))
    r=pd.DataFrame(rows);r.to_csv(out/'all_pathway_rank_effects.csv',index=False)
    for coll in ['Hallmark','Reactome']:
        r[r.collection.eq(coll)].pivot(index='pathway',columns='time_h',values='rank_biserial').to_csv(out/f'{coll.lower()}_timecourse_rank_effects.csv')
    pd.DataFrame([{'mapped_probes':len(merged),'genes':len(gene),'hallmark_sets':sum(x[0]=='Hallmark' for x in eligible),'reactome_sets':sum(x[0]=='Reactome' for x in eligible)}]).to_csv(out/'geo_qc_summary.csv',index=False)
    return len(gene),len(eligible)

def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,default=Path('data/external/GSE83654_series_matrix.txt.gz'));p.add_argument('--gpl',type=Path,default=Path('data/external/GPL13607-20416.txt'));p.add_argument('--hallmark',type=Path,default=Path('data/external/h.all.v2025.1.Hs.symbols.gmt'));p.add_argument('--reactome',type=Path,default=Path('data/external/c2.cp.reactome.v2025.1.Hs.symbols.txt'));p.add_argument('--out',type=Path,default=Path('outputs/transcriptomics'));p.add_argument('--download-geo',action='store_true');a=p.parse_args()
    if a.download_geo:download(URL_MATRIX,a.matrix)
    for path in [a.matrix,a.gpl,a.hallmark,a.reactome]:
        if not path.exists():raise FileNotFoundError(f'Missing external file {path}. See README for acquisition and licensing.')
    print('GEO result:',run(a.matrix,a.gpl,a.hallmark,a.reactome,a.out))
if __name__=='__main__':main()
