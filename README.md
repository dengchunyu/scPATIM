# scPATIM
Single-Cell Polygenic Analysis Pipeline for TIME Genetics

## Data Input: Get gwas data filterd by eqtls

```R
library(bigreadr)
library(readr)
library(dplyr)
GWAS_raw <-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/rawdata/finngen_R9_C3_DLBCL_EXALLC.gz")
colnames(GWAS_raw )<-c("chrom","pos","REF","ALT","rsid","nearest_genes","p","lp","beta","se","maf","maf_case","maf_control")
GWAS_raw$N<-287137
gwas_data<-GWAS_raw[!GWAS_raw$rsid %in% GWAS_raw$rsid[duplicated(GWAS_raw$rsid)],]
immsnp<-read_table('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_filter5e2.tsv',col_names = F)
snps<-intersect(gwas_data$rsid,immsnp$X2)
gwas_data<-gwas_data[gwas_data$rsid %in% snps,]
write.table(gwas_data,file= "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt",row.names=F,quote=F,sep="\t")
```

the input data of single cell has no difference with scPagwas.


## 2. Run the scPagwas

```bash
# before run the 1.Run_the_scPagwas.r, you should install the scPagwas and Seurat
source activate myr
cd /share/pub/dengcy/Cancer_Gwas/scPATIM/
scdata_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
gwas_ad='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt'
file_dir='/share/pub/dengcy/Cancer_Gwas/scPATIM/'
out='DLBCL_B'
Rscript 1.Run_the_scPagwas.r $scdata_ad $gwas_ad $file_dir $out
```

## 2. Get TRGs

```bash
Rscript 2.Get_TRGs.r $scdata_ad $file_dir $out
```

## 3. Get TRS and pvalue for each cell

Based on the scale of single-cell data and variations in memory computation, we can utilize either of the following two calculation methods.The results of these two methods are almost identical.

### 3.1 If you single cell is small.

```bash
Rscript 3.1.Get_TRS_pvalue_by_scPagwas.r $scdata_ad $file_dir $out 500 200
```

### 3.2 If you single cell is big.

```bash
#install the scDRS
source activate mypy
scRNA_h5ad_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad
gene_file=/share/pub/dengcy/Cancer_Gwas/scPATIM/DLBCL_B/_TRGs_PCC.csv
TRS_score_file=/share/pub/dengcy/Cancer_Gwas/scPATIM/DLBCL_B/DBCLC_B_scpatim.txt
celltype_file=/share/pub/dengcy/Cancer_Gwas/scPATIM/DLBCL_B/DBCLC_B_celltypes_pvalues.csv

python 3.2.Get_TRS_pvalue_by_scDRS.py --scRNA_h5ad_file $scRNA_h5ad_file --gene_file $gene_file --top_gene_num 500 --n_ctrl 200 --score_file $TRS_score_file --weight_pcc weight_pcc --group cell_type --celltype_file  $celltype_file
```

