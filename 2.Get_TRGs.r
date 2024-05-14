library(scPagwas)
library(Seurat)
Args <- commandArgs(T)

scdata = print(Args[1])
setwd = print(Args[2])
out = print(Args[3])

output.dirs<-out
output.prefix<-""

setwd(setwd)
load(paste0("./",out,"/Pagwas.RData"))
Single_data<-readRDS(scdata)

Random_PCC<- function(gPas,datamat,seed=1234,random_times=100,select_num=10000){
    message("Start to run the corSparse function in random!")
    set.seed(seed)
    sparse_cor_list<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(datamat),select_num)
    gPas_select <- data.matrix(gPas[index])
    sparse_cor <- scPagwas::corSparse(
      X = t(datamat[,index]),
      Y = gPas_select
    )
    colnames(sparse_cor) <- "PCC"
    sparse_cor[is.nan(sparse_cor)] <- 0
    sparse_cor[is.na(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
    }
    sparse_cor_list<-as.data.frame(sparse_cor_list)
    sparse_cor<- apply(sparse_cor_list,1,function(x) mean(x, na.rm = TRUE))
    return(data.matrix(sparse_cor))
}

scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
scPagwas.gPAS.score <- scPagwas::scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
Pagwas$data_mat <- Seurat::GetAssayData(Single_data, slot = "data", assay = "RNA")
Pagwas$data_mat <- scPagwas::as_matrix(Pagwas$data_mat)
pcc <- Random_PCC(gPas=scPagwas.gPAS.score,datamat=Pagwas$data_mat,seed=1234,random_times=100,select_num=floor(ncol(Pagwas$Pathway_single_results)/5))
colnames(pcc)<-'PCC'

# 提取前5%的scPagwas.gPAS.score
top5_idx <- order(scPagwas.gPAS.score, decreasing = TRUE)[1:(0.05 * length(scPagwas.gPAS.score))]

# 根据索引获取名称
top5_cellnames <- names(scPagwas.gPAS.score)[top5_idx]

# 在scPagwas.gPAS.score的名称中判断是否在top5_cellnames中
top5_idx2 <- names(scPagwas.gPAS.score) %in% top5_cellnames

# 对Pagwas$data_mat的每一行应用wilcox.test并进行bonferroni校正
p_list <- apply(Pagwas$data_mat, 1, function(x) {
    p.adjust(wilcox.test(x[top5_idx2], x[!top5_idx2], alternative = "greater")$p.value, method = "bonferroni")
})

# 创建pcc数据框并添加pvalue
pcc <- data.frame(pcc, pvalue = p_list)

# 对pvalue进行bonferroni校正
pcc <- within(pcc, {
    adj_pvalue <- p.adjust(pvalue, method = "bonferroni")
    adj_logp <- -log10(adj_pvalue)
    adj_logp <- ifelse(is.finite(-log10(adj_pvalue)), -log10(adj_pvalue), max(adj_logp, na.rm = TRUE) + 1)
})

pcc$weight_pcc<- pcc$adj_logp * pcc$PCC
utils::write.csv(pcc,file = paste0("./", output.dirs, "/",output.prefix,"_TRGs_PCC.csv"), quote = F)
