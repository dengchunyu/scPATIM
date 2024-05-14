library(scPagwas)
library(Seurat)
Args <- commandArgs(T)

scdata = print(Args[1])
#scdata<-'/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
filedir = print(Args[2])
#filedir<- '/share/pub/dengcy/Cancer_Gwas/scPATIM/'
out = print(Args[3])
#out <- 'DLBCL_B'
top_gene_num = print(Args[4])
top_gene_num <- as.numeric(top_gene_num)
n_ctrl = print(Args[5])
n_ctrl <- as.numeric(n_ctrl)


output.dirs<-out
output.prefix<-""

setwd(filedir)
load(paste0("./",out,"/Pagwas.RData"))
Single_data<-readRDS(scdata)

assay="RNA"

pcc <-read.csv(paste0("./", output.dirs, "/",output.prefix,"_TRGs_PCC.csv"),row.names=1)

scPagwas_topgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = T)[1:top_gene_num]]

Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes), name = c("scPagwas.TRS.Score"))

message("* Get Random Correct background pvalue for each single cell!")
correct_pdf <- scPagwas::Get_CorrectBg_p(Single_data=Single_data,
                                        scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                                        iters_singlecell=n_ctrl,
                                        n_topgenes=top_gene_num,
                                        scPagwas_topgenes=scPagwas_topgenes,
                                        assay=assay)

Pagwas$Random_Correct_BG_pdf <- correct_pdf

message("* Get Merged pvalue for each celltype!")
Merged_celltype_pvalue <- scPagwas::Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)

a <- data.frame(scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
Random_Correct_BG_p = correct_pdf$pooled_p,
Random_Correct_BG_adjp = correct_pdf$adj_p,
Random_Correct_BG_z = correct_pdf$pooled_z)
utils::write.csv(a,file = paste0("./", output.dirs, "/",output.prefix,"_singlecell_scPagwas_score_pvalue.Result.csv"),quote = F)

utils::write.csv(Merged_celltype_pvalue,
                file = paste0(
                    "./", output.dirs, "/", output.prefix,
                    "_Merged_celltype_pvalue.csv"
                ),
                quote = F
)



