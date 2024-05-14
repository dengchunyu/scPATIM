import scanpy as sc
import pandas as pd
import scdrs
import argparse
import numpy as np
import warnings
from scipy.stats import norm
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Get scDRS score')
parser.add_argument('--scRNA_h5ad_file', '-i', type=str, help='scRNA-seq h5ad file')
#scDRS gene and weight result file,default is scDRS_gene_weight.txt
parser.add_argument('--gene_file', '-g', type=str, help='scDRS gene and weight result file,default is gene_weight.txt')
#top gene number
parser.add_argument('--top_gene_num', '-t', type=int, help='top gene number')
#number of control
parser.add_argument('--n_ctrl', '-n', type=int, help='number of control')
#output:scDRS score file,default is scDRS_score.csv
parser.add_argument('--score_file', '-o', type=str, help='score file,default is scDRS_score.csv')
parser.add_argument('--weight_pcc', '-p', type=str, help='names for pcc or weight_pcc')
parser.add_argument('--group', '-u', type=str, help='the names for celltype group')
parser.add_argument('--celltype_file', '-f', type=str, help='result file for celltype result')

args = parser.parse_args()

adata_file = args.scRNA_h5ad_file
scDRS_gene_file = args.gene_file
top_gene_num = args.top_gene_num
n_ctrl = args.n_ctrl
scDRS_score_file = args.score_file
weight_col = args.weight_pcc
group = args.group
celltype_file = args.celltype_file

adata = sc.read_h5ad(adata_file)
scdrs.preprocess(adata)
if scDRS_gene_file.endswith('.txt'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep='\t')
elif scDRS_gene_file.endswith('.csv'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep=',')

#判断scDRS_gene_file是否有gene，如果没有判断是否有gene_symbol，如果有gene_symbol则将其改为gene，如果没有则报错
if 'GENE' not in scDRS_gene.columns:
    if 'gene_symbol' in scDRS_gene.columns:
        scDRS_gene = scDRS_gene.rename(columns={'gene_symbol':'GENE'})
    else:
        scDRS_gene = scDRS_gene.rename(columns={'Unnamed: 0':'GENE'})
    

scDRS_gene = scDRS_gene.rename(columns={weight_col:'weight'})
scDRS_gene = scDRS_gene.dropna(subset=['weight'])
max_weight  = scDRS_gene['weight'].replace([np.inf, -np.inf], np.nan).max()
scDRS_gene['weight'] = scDRS_gene['weight'].replace([np.inf, -np.inf], max_weight + 1)
scDRS_gene = scDRS_gene.sort_values(by='weight',ascending=False).iloc[0:top_gene_num,:]

#计算scDRS得分
df_score = scdrs.score_cell(data=adata,gene_list=scDRS_gene['GENE'],gene_weight=scDRS_gene['weight'],ctrl_match_key="mean_var",n_ctrl=n_ctrl,weight_opt="vs",return_ctrl_raw_score=False,return_ctrl_norm_score=True,verbose=False)
df_score.to_csv(scDRS_score_file)

def p_merge(pvalues):
    zvalues = -np.sqrt(2) * norm.ppf(pvalues / 2)
    ztotal = np.mean(zvalues)
    p_total = norm.cdf(-abs(ztotal))
    return p_total

def merge_celltype_p(single_p, celltype):
    celltype_p = pd.DataFrame({'celltype': celltype, 'pvalue': single_p})
    celltype_p = celltype_p.groupby('celltype')['pvalue'].agg(p_merge).reset_index()
    return celltype_p
if not args.celltype_file:
    pass
else:
    drs_df=df_score.iloc[:,[1,3]]
    drs_df.columns='scDRS_'+drs_df.columns
    drs_df['Significant_cells'] = drs_df['scDRS_pval'].apply(lambda x: 1 if x < 0.05 else 0)
    drs_df.index=adata.obs.index
    merged_df = pd.concat([adata.obs, drs_df], axis=1)
    adata.obs=merged_df
    adata.write(adata_file)
    ct_value=merge_celltype_p(merged_df['scDRS_pval'], merged_df[group])
    ct_value.to_csv(celltype_file)
