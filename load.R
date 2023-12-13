suppressPackageStartupMessages({
    library(readr)
    library(stringr)
    library(dplyr)
    library(magrittr)
    library(SingleCellExperiment)
    library(scater)
    library(flexmix)
    library(splines)
    library(BiocParallel)
    library(biomaRt)
    library(miQC)
    library(Seurat)
    library(SeuratDisk)
})
runs <-
    list.dirs() %>%
    .[stringr::str_detect(string = ., pattern = "GSM440413")] %>%
    str_remove(pattern = "./")

prepRun <- function(pathr) {
    renm <- c("P30_VPH_M2", "P30_VPH_F2", "P30_VPH_F3", "P30_VPH_M3")
    names(renm) <- c("GSM4404135_AJ18003", "GSM4404136_AJ18004",
                     "GSM4404137_AJ19001", "GSM4404138_AJ19002")
    prj <- pathr
    if (prj %in% names(renm)) {
        prj <- renm[prj]
        sex <- if_else(str_detect(prj, "M"), "M", "F")
        tech <- if_else(str_detect(prj, "3$"), "10xv3", "10xv2")
    }
    result_f <- sprintf("%s.h5Seurat", prj)
    mtx <- Read10X(data.dir = sprintf("%s/", pathr))
    ## Initialize the Seurat object with the raw (non-normalized data).
    srt <- CreateSeuratObject(
        counts = mtx,
        project = prj,
        min.cells = 0,
        min.features = 200
    )

    srt@meta.data %<>%
        tibble::rownames_to_column(var = "bc_name")
    rownames(srt@meta.data) <- colnames(srt)
    srt$age <- "P30"
    srt$sex <- sex
    srt$tech <- tech
    srt$study_id <- "mickelsen_2020"

    return(srt)
}
nameRun <- function(pathr) {
    renm <- c("P30_VPH_M2", "P30_VPH_F2", "P30_VPH_F3", "P30_VPH_M3")
    names(renm) <- c("GSM4404135_AJ18003", "GSM4404136_AJ18004",
                     "GSM4404137_AJ19001", "GSM4404138_AJ19002")
    prj <- pathr
    if (prj %in% names(renm)) {
        prj <- renm[prj]
    }
    return(prj)
}
srt_list <- runs %>% purrr::map(prepRun)
names(srt_list) <- runs %>% purrr::map_chr(nameRun)
## 1
sce <- as.SingleCellExperiment(srt_list[[1]])
mt_genes <- grepl("^mt-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
feature_ctrls
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

plotMetrics(sce)
plotColData(sce, x = "sum", y="detected", colour_by="subsets_mito_percent")
plotColData(sce, x = "sum", y="subsets_mito_percent", colour_by="detected")
sce <- sce[, sce$subsets_mito_percent < 15]
sce <- sce[, sce$detected > 1500]
sce <- sce[, sce$detected < 6500]
sce <- sce[, sce$sum < 20000]
srt_list[[1]] <- as.Seurat(sce)
rm(sce, model, model2, feature_ctrls, mt_genes)


## 2
sce <- as.SingleCellExperiment(srt_list[[2]])
mt_genes <- grepl("^mt-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
feature_ctrls
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

plotMetrics(sce)
model <- mixtureModel(sce)
plotModel(sce, model)
plotFiltering(sce, model)
plotColData(sce, x = "sum", y="detected", colour_by="subsets_mito_percent")
plotColData(sce, x = "sum", y="subsets_mito_percent", colour_by="detected")
sce <- sce[, sce$subsets_mito_percent < 15]
sce <- sce[, sce$detected > 1500]
sce <- sce[, sce$detected < 7000]
sce <- sce[, sce$sum < 20000]
srt_list[[2]] <- as.Seurat(sce)
rm(sce, model, model2, feature_ctrls, mt_genes)


## 3
sce <- as.SingleCellExperiment(srt_list[[3]])
mt_genes <- grepl("^mt-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
feature_ctrls
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

plotMetrics(sce)
plotColData(sce, x = "sum", y="detected", colour_by="subsets_mito_percent")
plotColData(sce, x = "sum", y="subsets_mito_percent", colour_by="detected")
sce <- sce[, sce$subsets_mito_percent < 15]
sce <- sce[, sce$detected > 2000]
sce <- sce[, sce$detected < 7500]
sce <- sce[, sce$sum < 40000]
srt_list[[3]] <- as.Seurat(sce)
rm(sce, model, model2, feature_ctrls, mt_genes)


## 4
sce <- as.SingleCellExperiment(srt_list[[4]])
mt_genes <- grepl("^mt-",  rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
feature_ctrls
sce <- addPerCellQC(sce, subsets = feature_ctrls,
                    BPPARAM = BiocParallel::MulticoreParam())

plotMetrics(sce)
plotColData(sce, x = "sum", y="detected", colour_by="subsets_mito_percent")
plotColData(sce, x = "sum", y="subsets_mito_percent", colour_by="detected")
sce <- sce[, sce$subsets_mito_percent < 15]
sce <- sce[, sce$detected > 2000]
sce <- sce[, sce$detected < 7500]
sce <- sce[, sce$sum < 40000]
srt_list[[4]] <- as.Seurat(sce)
rm(sce, model, model2, feature_ctrls, mt_genes)

mickelsen2020_combined_vph <-
    merge(srt_list[["P30_VPH_M2"]],
          y = c(srt_list[["P30_VPH_F2"]],
                srt_list[["P30_VPH_F3"]],
                srt_list[["P30_VPH_M3"]]),
          add.cell.ids = c("P30_VPH_M2", "P30_VPH_F2", "P30_VPH_F3", "P30_VPH_M3"),
          project = "PVH")
glimpse(mickelsen2020_combined_vph@meta.data)
table(mickelsen2020_combined_vph$orig.ident)
SaveH5Seurat(mickelsen2020_combined_vph, filename = "mickelsen2020_vph.h5Seurat")
Convert("mickelsen2020_vph.h5Seurat", dest = "h5ad")


