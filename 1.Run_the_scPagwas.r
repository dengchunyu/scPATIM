if ("scPagwas" %in% rownames(installed.packages())) {
    package_version <- packageVersion("scPagwas")
    if (package_version < "1.3") {
        stop("Error: scPagwas version must be greater than 1.3.")
    }
} else {
    stop("Error: scPagwas package is not installed.")
}


library(Seurat)
library(scPagwas)
Args <- commandArgs(T)

scdata = print(Args[1])
gwas = print(Args[2])
file_dir = print(Args[3])
out = print(Args[4])

setwd(file_dir)
if (!file.exists("mp_immune_genelst.RData")) {
    stop("Error: mp_immune_genelst.RData file not found in the current working directory.")
}
load("mp_immune_genelst.RData")
output.dirs<-out
output.prefix<-""

Pagwas<-scPagwas::scPagwas_main(Pagwas =NULL,
        Single_data =scdata,
        gwas_data =gwas,
        output.prefix="",
        output.dirs=out,
        n_topgenes = 500,
        Pathway_list=mp_immune_genelst,
        iters_singlecell = 0,
        seurat_return=F,
        celltype = FALSE,
        block_annotation = block_annotation,
        chrom_ld = chrom_ld)

save(Pagwas,file=paste0("./",out,"/Pagwas.RData"))