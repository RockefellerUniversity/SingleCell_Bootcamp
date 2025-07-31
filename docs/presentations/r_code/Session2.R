params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))

suppressPackageStartupMessages(require(SeuratWrappers))
suppressPackageStartupMessages(require(ggplot2))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----warning=F, message=F, echo=F, eval=F-------------------------------------
# library(RCurl)
# if(!url.exists("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz")){stop("Download path broken")}
# 


## ----eval=F-------------------------------------------------------------------
# download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz","pbmc8k_filtered_gene_bc_matrices.tar.gz")


## ----eval=F-------------------------------------------------------------------
# untar("pbmc8k_filtered_gene_bc_matrices.tar.gz")
# 


## ----eval=F-------------------------------------------------------------------
# dir()
# 


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Ranger to Seurat

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger to Seurat


---
"    
  )
  
}



## ----echo=T,eval=TRUE---------------------------------------------------------
mtx_dir <- "filtered_gene_bc_matrices/GRCh38"




## ----load_data,include=TRUE,eval=F--------------------------------------------
# library(Seurat)
# 
# mtx <- Seurat::Read10X(mtx_dir)
# 


## ----eval=F, echo=F-----------------------------------------------------------
# a <- is(mtx)
# b <- head(mtx)
# save(a,b, file="data/pbmc8k_mtx_res.RData")


## ----eval=F, echo=F-----------------------------------------------------------
# library(Seurat)
# save(mtx, file = "data/pbmc8k_seurat_read.RData")
# load("data/pbmc8k_seurat_read.RData")


## ----eval=T, echo=F-----------------------------------------------------------
library(Seurat)
load("data/pbmc8k_mtx_res.RData")


## ----eval=F, echo=T-----------------------------------------------------------
# is(mtx)


## ----eval=T, echo=F-----------------------------------------------------------
a


## ----eval=F, echo=T-----------------------------------------------------------
# head(mtx)


## ----eval=T, echo=F-----------------------------------------------------------
b


## ----load_h5,include=TRUE,eval=FALSE------------------------------------------
# h5_file <- "path to matrix h5 file"
# h5_file <- "~/Downloads/pbmc8k_raw_gene_bc_matrices_h5.h5"
# mtx <- Seurat::Read10X_h5(h5_file)


## ----load_CreateOBJ,include=TRUE,eval=TRUE------------------------------------
sample_id <- "PBMC_8k" # sample name
min_gene <- 200
min_cell <- 10 


## ----eval=F-------------------------------------------------------------------
# seu_obj <- Seurat::CreateSeuratObject(mtx, project=sample_id, min.cells=min_cell, min.features=min_gene)


## ----eval=F, echo=F-----------------------------------------------------------
# save(seu_obj, file="data/pbmc8k_seu_obj_raw.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/pbmc8k_seu_obj_raw.RData") 


## -----------------------------------------------------------------------------
seu_obj[["dset"]] <- sample_id # Create a category for sample
seu_obj <- Seurat::RenameCells(seu_obj, add.cell.id=sample_id) # add sample name in front of cell barcode


## ----laod_CreatOBJ_pres,include=TRUE,eval=TRUE--------------------------------
seu_obj


## -----------------------------------------------------------------------------
head(seu_obj,2)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Normalization, feature selection, and data scaling

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Normalization, feature selection, and data scaling


---
"    
  )
  
}



## ----norm_log,include=TRUE,eval=TRUE------------------------------------------
seu_obj <- NormalizeData(seu_obj, normalization.method="LogNormalize")


## -----------------------------------------------------------------------------
seu_obj <- FindVariableFeatures(seu_obj, select.method="vst", nfeatures=3000)


## -----------------------------------------------------------------------------
seu_obj <- ScaleData(seu_obj)


## ----norm_plotHVF,include=TRUE, eval=F----------------------------------------
# top10 <- head(VariableFeatures(seu_obj), 10)
# 
# plot1 <- VariableFeaturePlot(seu_obj)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2


## ----echo=F,fig.height=4,fig.width=7------------------------------------------
top10 <- head(VariableFeatures(seu_obj), 10)

plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


## ----eval=F, echo=F-----------------------------------------------------------
# download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz","pbmc8k_filtered_gene_bc_matrices.tar.gz")
# untar("pbmc8k_filtered_gene_bc_matrices.tar.gz")


## ----include=F,echo=T,eval=F--------------------------------------------------
# library(Seurat)
# 
# sample_id <- "PBMC_8k"
# min_gene <- 200
# min_cell <- 10
# 
# mtx_dir <- "filtered_gene_bc_matrices/GRCh38"
# seu_obj <- Seurat::Read10X(mtx_dir)
# seu_obj <- Seurat::CreateSeuratObject(seu_obj, project=sample_id, min.cells=min_cell, min.features=min_gene)
# seu_obj <- NormalizeData(seu_obj, normalization.method="LogNormalize")
# seu_obj <- FindVariableFeatures(seu_obj, select.method="vst", nfeatures=3000)
# seu_obj <- ScaleData(seu_obj)


## ----eval=F, echo=F-----------------------------------------------------------
# unlink("pbmc8k_filtered_gene_bc_matrices.tar.gz", recursive=TRUE)
# unlink("filtered_feature_bc_matrix/GRCh38", recursive=TRUE)


## ----norm_sct, include=TRUE, eval=F-------------------------------------------
# seu_obj <- SCTransform(seu_obj, variable.features.n = 3000)
# seu_obj


## ----eval=F, echo=F-----------------------------------------------------------
# SCT_assay <- seu_obj@assays$SCT
# save(SCT_assay, file="data/pbmc8k_SCT_assay.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/pbmc8k_SCT_assay.RData")
seu_obj@assays$SCT <- SCT_assay
DefaultAssay(seu_obj) <- "SCT"
seu_obj


## ----include=TRUE, eval=T-----------------------------------------------------
DefaultAssay(seu_obj) <- "RNA"

DefaultAssay(seu_obj) <- "SCT"


## ----include=F,warning=F, message=F, echo=FALSE-------------------------------
rm(SCT_assay)
gc()



## ----norm_comp,include=TRUE,eval=F--------------------------------------------
# log_mat <- GetAssayData(seu_obj,assay="RNA",slot="data")
# log_mat <- as.matrix(log_mat)
# log_avgExp <- rowMeans(log_mat)
# 
# sct_mat <- GetAssayData(seu_obj,assay="SCT",slot="data")
# sct_mat <- as.matrix(sct_mat)
# sct_avgExp <- rowMeans(sct_mat)


## ----eval=F, echo=F-----------------------------------------------------------
# 
# exp_sct <- list(log_avgExp,sct_avgExp)
# names(exp_sct) <- c("logNorm","SCTransform")
# save(exp_sct, file="data/pbmc8k_exp_sct.RData")
# 


## ----eval=T, echo=F-----------------------------------------------------------

load("data/pbmc8k_exp_sct.RData")
log_avgExp <- exp_sct$logNorm
sct_avgExp <- exp_sct$SCTransform



## ----eval=F-------------------------------------------------------------------
# library(ggplot2)
# 
# dat <- data.frame(logNorm=log_avgExp, SCT=sct_avgExp)
# cor_val <- cor.test(log_avgExp,sct_avgExp,method = "spearman")
# 
# ggplot(dat,aes(x=logNorm,y=SCT))+geom_point()+geom_smooth()+
#   labs(x="Log_Normalization",y="SCTransform",subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
#   theme_classic()


## ----echo=F, fig.height=4,fig.width=7-----------------------------------------
dat <- data.frame(logNorm=log_avgExp, SCT=sct_avgExp)
cor_val <- cor.test(log_avgExp,sct_avgExp,method = "spearman")

ggplot(dat,aes(x=logNorm,y=SCT))+geom_point()+geom_smooth()+
  labs(x="Log_Normalization",y="SCTransform",subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
  theme_classic()


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(dat, log_mat, sct_mat, plot1, plot2, SCT_assay)
DefaultAssay(seu_obj) <- "RNA"
gc()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Dimension Reduction

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html>

---
"
  )
}else{
  cat("# Dimension Reduction


---
"
  )

}



## ----dimRed_pca,include=TRUE,eval=T-------------------------------------------
set.seed(1001)
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- RunPCA(seu_obj, assay = "RNA", npcs = 50)



## ----dimRed_elbow,include=TRUE,eval=TRUE, fig.height=4,fig.width=7------------
ElbowPlot(seu_obj, ndims = 50, reduction = "pca")


## ----dimRed_elbow2,include=TRUE,eval=T----------------------------------------
library(dplyr)
# Determine percent of variation associated with each PC
pct <- seu_obj[["pca"]]@stdev / sum(seu_obj[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)


## ----include=TRUE,eval=T------------------------------------------------------
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC, last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pc <- min(co1, co2)

pc



## ----dimRed_elbow2_rep,include=TRUE,eval=F------------------------------------
# plot_df <- data.frame(pct = pct, cumu = cumu,rank = 1:length(pct))
# 
# ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
#   geom_vline(xintercept = 90, color = "grey") +
#   geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#   theme_bw()


## ----eval=T, echo=F, fig.height=4,fig.width=7---------------------------------
plot_df <- data.frame(pct = pct, cumu = cumu,rank = 1:length(pct))

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


## ----dimRed_UMAP,include=TRUE,eval=T------------------------------------------
seu_obj <- FindNeighbors(seu_obj,dims = 1:pc, reduction = "pca")

seu_obj <- RunUMAP(seu_obj,dims = 1:pc, reduction = "pca")


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_obj)



## ----include=TRUE,eval=T------------------------------------------------------
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
seu_obj[["cluster_byDefault"]] <- seu_obj$seurat_clusters


## ----include=TRUE,eval=T------------------------------------------------------
DimPlot(seu_obj, group.by = "cluster_byDefault")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Droplet processing

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Droplet processing

---
"    
  )
  
}



## input_h5=the_raw_matrix_in_h5_format_from_cellranger #essential
## output_h5=assign_the_h5_file_path_for_the_cellbender_corrected_matrix # essential
## fpr=threshold_of_FALSE_POSITIVE_RATE # default 0.01
## epochs=number_of_epochs_to_train # default 150
## num_train=number_of_times_to_attempt_to_train_the_model # default 1. would speed up while setting greater
## 
## cellbender remove-background --input $input_h5 --output $output_h5 --fpr $fpr --epochs $epochs --num-training-tries $num_train --cuda

## ptrepack --complevel 5 celbender_filtered.h5:/matrix celbender_filtered_forseurat.h5:/matrix
## 

## ----sec3_dropProc_cbFilt_loadData,include=F,eval=T, echo=FALSE---------------


load("data/pbmc8k_cbFilt_20250131_filtered_mtx.RData") 



## ----include = T, eval = F, echo = T------------------------------------------
# 
# 
# cbFilt_mtx <- Read10X_h5("~/Downloads/cbFilt_PBMC8K_20250131_filtered_forseurat.h5")
# 


## ----sec_mergeData_funcUsed_dataProc,include=TRUE-----------------------------
data_proc <- function(seu){
  seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, select.method="vst", nfeatures=2000)
  return(seu)}


## ----sec3_mergeData_funcUsed_quickClust,include=TRUE--------------------------
quick_clust <- function(seu){
  set.seed(42)
  seu <- RunPCA(seu, npcs=30, verbose=FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, verbose=FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:10, verbose=FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose=FALSE)
  return(seu)}


## ----sec3_dropProc_cb_Filt_dataPRoc,include=TRUE,eval=T, message=F, warning=F----

message("processing matrix from CellBender")
seu <- CreateSeuratObject(cbFilt_mtx)
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
seu_cbFilt <- seu


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu)
gc()
  


## ----fig.height=4,fig.width=6-------------------------------------------------
DimPlot(seu_obj ,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## ----fig.height=4,fig.width=6-------------------------------------------------
DimPlot(seu_cbFilt,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## ----fig.height=4,fig.width=7, warning=F--------------------------------------

mark_gene <- c("LYZ","HLA-DRA")

FeaturePlot(seu_obj,features = mark_gene,pt.size = 0)



## ----fig.height=3,fig.width=7, warning=F--------------------------------------


VlnPlot(seu_obj,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)



## ----fig.height=4,fig.width=7, warning=F--------------------------------------

FeaturePlot(seu_cbFilt,features = mark_gene,pt.size = 0)



## ----fig.height=3,fig.width=7, warning = F------------------------------------

VlnPlot(seu_cbFilt,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## ----echo=F, warning=F, message=FALSE-----------------------------------------
rm(cbFilt_mtx, filt_mtx, seu_obj)
gc()


## -----------------------------------------------------------------------------
sample_id <- "PBMC_8k"
seu_cbFilt[["dset"]] <- sample_id # Create a category for sample
seu_cbFilt <- Seurat::RenameCells(seu_cbFilt, add.cell.id=sample_id) # add sample name in front of cell barcode



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Mitochondrial Proportion

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Mitochondrial Proportion


---
"    
  )
  
}



## -----------------------------------------------------------------------------
mt_gene <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
             "MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8","MT-CYB")
mt_gene_det <- mt_gene[mt_gene %in% rownames(seu_cbFilt)]
seu_cbFilt[["percent.mt"]] <- PercentageFeatureSet(seu_cbFilt, features = mt_gene_det)
summary(seu_cbFilt$percent.mt)


## ----load_estMT,include=TRUE,eval=TRUE----------------------------------------
seu_cbFilt[["percent.mt2"]] <- PercentageFeatureSet(seu_cbFilt, pattern = "^MT-")
summary(seu_cbFilt$percent.mt2)



## ----eval =F------------------------------------------------------------------
# library(ggplot2)
# dat <- data.frame(byPattern=seu_cbFilt$percent.mt, byGene=seu_cbFilt$percent.mt2,stringsAsFactors = FALSE)
# cor_val <- cor.test(dat$byPattern,dat$byGene,method = "spearman")
# ggplot(dat,aes(x=byPattern,y=byGene))+geom_point()+geom_smooth()+
#   labs(x="% of MT, est by genes",y="% of MT,est by pattern",
#        subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
#   theme_classic()


## ----echo =F, fig.height=4,fig.width=7----------------------------------------
library(ggplot2)
dat <- data.frame(byPattern=seu_cbFilt$percent.mt, byGene=seu_cbFilt$percent.mt2,stringsAsFactors = FALSE)
cor_val <- cor.test(dat$byPattern,dat$byGene,method = "spearman")
ggplot(dat,aes(x=byPattern,y=byGene))+geom_point()+geom_smooth()+
  labs(x="% of MT, est by genes",y="% of MT,est by pattern",
       subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
  theme_classic()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Cycle Phases

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Cycle Phases


---
"    
  )
  
}



## ----ccPhase_Seurat,include=TRUE,eval=TRUE------------------------------------
feat_s <- cc.genes$s.genes
feat_g2m <- cc.genes$g2m.genes

feat_s
feat_g2m


## ----ccPhase_plot_Seurat,include=TRUE,eval=TRUE-------------------------------
DefaultAssay(seu_cbFilt) <- "RNA"
seu_cbFilt <- CellCycleScoring(seu_cbFilt, s.features = feat_s, g2m.features = feat_g2m)



## -----------------------------------------------------------------------------

dat_s <- data.frame(cell_id=Cells(seu_cbFilt), cat="S_Score", Phase=seu_cbFilt$Phase, score=seu_cbFilt$S.Score)

dat_g2m <- data.frame(cell_id=Cells(seu_cbFilt), cat="G2M_Score", Phase=seu_cbFilt$Phase, score=seu_cbFilt$G2M.Score)

dat <- rbind(dat_s, dat_g2m)
dat$Phase <- factor(dat$Phase, levels = c("G1","S","G2M"))


## ----fig.height=4,fig.width=7-------------------------------------------------

ggplot(dat,aes(x=Phase, y=score, fill=Phase))+geom_violin()+
  labs(x="",y="Score",fill="Phase")+
  facet_wrap(~cat)+theme_classic()


## ----ccPhase_cyclon_prep,include=TRUE,eval=TRUE, warning=F, message=F---------
library(scran)

sce <- as.SingleCellExperiment(seu_cbFilt, assay = "RNA")
rowData(sce)$SYMBOL <- rownames(sce)


## -----------------------------------------------------------------------------
sce


## -----------------------------------------------------------------------------

load("data/ccGene_mouse_human_geneSymbol_ensemblID_20220817.RData")
ccGene_hs <- ccGene_mm_hs$human_symbol
lapply(ccGene_hs, function(x){head(x,2)})


## ----ccPhase_cyclon_proc,include=TRUE,eval=FALSE------------------------------
# assignments <- cyclone(sce, ccGene_hs, gene.names=rowData(sce)$SYMBOL)
# 


## ----echo=F, eval=F-----------------------------------------------------------
# save(assignments, file="data/pbmc8k_cyclone.RData")


## ----echo=F-------------------------------------------------------------------

load("data/pbmc8k_cyclone.RData")


## -----------------------------------------------------------------------------
lapply(assignments, head)


## -----------------------------------------------------------------------------
seu_cbFilt[["cyclon_Phase"]] <- assignments$phases
seu_cbFilt[["cyclon_G1Score"]] <- assignments$scores$G1
seu_cbFilt[["cyclon_SScore"]] <- assignments$scores$S
seu_cbFilt[["cyclon_G2MScore"]] <- assignments$scores$G2M


## ----ccPhase_cyclon_boxPlot,include=TRUE,eval=T-------------------------------
dat_g1 <- data.frame(cell_id=Cells(seu_cbFilt), cat="cyclon_G1Score", Phase=seu_cbFilt$cyclon_Phase, score=seu_cbFilt$cyclon_G1Score)

dat_s <- data.frame(cell_id=Cells(seu_cbFilt), cat="cyclon_SScore", Phase=seu_cbFilt$cyclon_Phase, score=seu_cbFilt$cyclon_SScore)

dat_g2m <- data.frame(cell_id=Cells(seu_cbFilt), cat="cyclon_G2MScore", Phase=seu_cbFilt$cyclon_Phase, score=seu_cbFilt$cyclon_G2MScore)

dat <- rbind(dat_g1,dat_s,dat_g2m)

dat$Phase <- factor(dat$Phase,levels = c("G1","S","G2M"))
dat$cat <- factor(dat$cat,levels = c("cyclon_G1Score","cyclon_SScore","cyclon_G2MScore"))


## ----fig.height=4,fig.width=7, warning=F--------------------------------------
ggplot(dat,aes(x=Phase,y=score,fill=Phase))+geom_violin()+
  labs(x="",y="Score",fill="Phase")+
  facet_wrap(~cat)+theme_classic()


## -----------------------------------------------------------------------------
table(seu_cbFilt$Phase)


## -----------------------------------------------------------------------------
table(seu_cbFilt$cyclon_Phase)


## -----------------------------------------------------------------------------
table(seu_cbFilt$Phase,seu_cbFilt$cyclon_Phase)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Detecting Doublets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Detecting Doublets


---
"    
  )
  
}



## ----eval=F, echo=FALSE-------------------------------------------------------
#  library(Herper)
# 
#  conda_install  <- install_CondaTools("scrublet", "scRNA", pathToMiniConda = "../mini")
# 
#  Sys.setenv('RETICULATE_PYTHON'=file.path(conda_install$pathToEnvBin, "python"))


## ----eval=F-------------------------------------------------------------------
# library(reticulate)
# reticulate::py_install("scrublet")
# 
# scr <- reticulate::import("scrublet")
# 


## ----det_doublet_est,include=TRUE,eval=F--------------------------------------
# 
# mat <- GetAssayData(seu_cbFilt, assay = "RNA", slot = "counts")
# mat <- as.matrix(mat)
# 


## ----eval=F-------------------------------------------------------------------
# # Loading the scrublet library
# scr <- import("scrublet")
# # Run scrublet
# scrub <- scr$Scrublet(t(mat))
# # Extract scrublet results
# doublet <- scrub$scrub_doublets()
# names(doublet) <- c("doublet_score","doublet")


## ----eval=F, echo=F-----------------------------------------------------------
# save(doublet,file = "data/pbmc8k_doublet.RData")
# 


## ----eval=T, echo=F-----------------------------------------------------------

load("data/pbmc8k_doublet.RData")



## -----------------------------------------------------------------------------
summary(doublet$doublet_score)


## -----------------------------------------------------------------------------
table(doublet$doublet)


## ----det_doublet_pres,include=TRUE,eval=TRUE----------------------------------
seu_cbFilt[["doublet_score"]] <- doublet$doublet_score
seu_cbFilt[["doublet"]] <- doublet$doublet


## ----fig.height=4,fig.width=7-------------------------------------------------

VlnPlot(seu_cbFilt, group.by = "doublet",
        features = c("doublet_score", "nCount_RNA","nFeature_RNA"),
        pt.size = 0)


## ----fig.height=4,fig.width=7-------------------------------------------------
FeatureScatter(seu_cbFilt,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",pt.size = 0.1,group.by = "doublet")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quality Assessment

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quality Assessment


---
"    
  )
  
}



## ----qcPlot_vlnPlot_pres1,include=TRUE,eval=F---------------------------------
# VlnPlot(seu_cbFilt, group.by = "dset",
#         features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#         pt.size = 0)


## ----include=TRUE,echo=F, fig.height=4,fig.width=7----------------------------
VlnPlot(seu_cbFilt, group.by = "dset",
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        pt.size = 0)


## ----include=TRUE,eval=T, echo=F, fig.height=3,fig.width=6--------------------

FeatureScatter(seu_cbFilt,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",
               group.by = "doublet")


## ----include=TRUE,eval=F------------------------------------------------------
# FeatureScatter(seu_cbFilt, feature1 = "nCount_RNA", feature2 = "percent.mt")


## ----qcPlot_scatter_pres2,include=TRUE,eval=TRUE, echo=F, fig.height=4,fig.width=7----
FeatureScatter(seu_cbFilt, feature1 = "nCount_RNA", feature2 = "percent.mt")


## ----include=TRUE,eval=F------------------------------------------------------
# RidgePlot(seu_cbFilt,group.by = "doublet",features = c("doublet_score"))


## ----qcPlot_ridgePlot_pres,include=TRUE,eval=TRUE, echo=F, fig.height=3,fig.width=6----
RidgePlot(seu_cbFilt,group.by = "doublet",features = c("doublet_score"))


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Filtering debris and doublets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Filtering debris and doublets


---
"    
  )
  
}



## ----filtCell_pres,include=TRUE,eval=TRUE-------------------------------------
table(seu_cbFilt$doublet=="TRUE" | seu_cbFilt$percent.mt >= 10)



## ----eval=T-------------------------------------------------------------------
seu_filt <- subset(seu_cbFilt, subset=doublet=="FALSE" & 
                     percent.mt < 10)



## ----eval = T, echo = F-------------------------------------------------------

VlnPlot(seu_cbFilt, group.by = "dset",
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        pt.size = 0)




## ----eval = T, echo = F-------------------------------------------------------

VlnPlot(seu_filt, group.by = "dset",
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        pt.size = 0)




## ----echo=F, eval=TRUE--------------------------------------------------------
#save(seu_filt,file="data/pbmc8k_seu_filt.RData")
rm(doublet)
rm(seu_cbFilt)
#load("data/pbmc8k_seu_filt.RData")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Regress out confounders

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Regress out confounders


---
"    
  )
  
}



## ----filtCell_scaleData,include=TRUE,eval=F-----------------------------------
# 
# #First let's replace NA values in the cyclon phase columns because ScaleData doesn't like this
# cyclone_cols_numeric <- c("cyclon_G1Score", "cyclon_G2MScore", "cyclon_SScore")
# # For each numeric column, replace NA with 0
# for (col_name in cyclone_cols_numeric) {
#  seu_filt@meta.data[[col_name]][is.na(seu_filt@meta.data[[col_name]])] <- 0
# }
# 
# pot_conf <- c("percent.mt","doublet_score","cyclon_G1Score","cyclon_SScore","cyclon_G2MScore")
# seu_filt <- ScaleData(seu_filt, vars.to.regress = pot_conf)
# 


## ----eval = T, warning=F, message=F-------------------------------------------

seu_filt <- data_proc(seu_filt)



## ----eval = T, warning=F, message=F-------------------------------------------

seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)
# Determine percent of variation associated with each PC
pct <- seu_filt[["pca"]]@stdev / sum(seu_filt[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC, last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pc <- min(co1, co2)



## ----eval = T, warning=F, message=F-------------------------------------------

seu_filt <- FindNeighbors(seu_filt,dims = 1:pc, reduction = "pca")
seu_filt <- RunUMAP(seu_filt,dims = 1:pc, reduction = "pca")



## ----eval=F, echo=F-----------------------------------------------------------
# 
# save(seu_filt, file="data/pbmc8k_SCT3.RData")
# saveRDS(seu_filt, file="data/pbmc8k_SCT3.rds")
# save(sce, file="data/pbmc8k_sce_updated.RData")
# 


## ----sceload, eval=T, echo=F--------------------------------------------------
load("data/pbmc8k_sce_updated.RData")


## ----sctload, eval=T, echo=F--------------------------------------------------
load("data/pbmc8k_SCT3.RData")
#seu_filt <- readRDS("data/pbmc8k_SCT3.rds")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Clustering

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html>

---
"
  )
}else{
  cat("# Clustering


---
"
  )
  
}



## ----include=TRUE,eval=T------------------------------------------------------
seu_filt <- FindClusters(seu_filt, resolution = 0.5)
seu_filt[["cluster_byDefault"]] <- seu_filt$seurat_clusters


## ----fig.height=4,fig.width=7-------------------------------------------------

DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2)+NoLegend()




## ----clust_resoEval,include=TRUE,eval=T, warning=FALSE, message=FALSE---------
library(clustree)

reso <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
reso_res <- lapply(1:length(reso), function(x,seu_filt,reso){
  seu_filt <- FindClusters(seu_filt,resolution = reso[x])
  clust <- setNames(seu_filt$seurat_clusters,Cells(seu_filt))
  return(clust)}, seu_filt, reso)
names(reso_res) <- paste0("k",1:length(reso))


## ----eval=F-------------------------------------------------------------------
# remotes::install_github("thomasp85/tweenr")


## ----include=TRUE,eval=T------------------------------------------------------
k_tab <- do.call(cbind,reso_res)
k_dat <- as.data.frame(k_tab)

head(k_dat,2)


## ----clust_resoEval2,include=TRUE,eval=F--------------------------------------
# clustree(k_dat, prefix = "k", node_colour = "sc3_stability")


## ----include=TRUE,eval=T, echo=F, warning=FALSE, fig.height=4,fig.width=7-----
clustree(k_dat, prefix = "k", node_colour = "sc3_stability")


## ----clust_resoFix3,include=TRUE,eval=T, fig.height=4,fig.width=7, warning=FALSE, message=FALSE----
seu_filt <- FindClusters(seu_filt, resolution = 0.3)


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Optimized Clusters")


## ----include=TRUE,eval=T, fig.height=4,fig.width=7----------------------------
DimPlot(seu_filt, group.by = "cluster_byDefault",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Default Clusters")


## ----fig.height=4,fig.width=6-------------------------------------------------

FeaturePlot(seu_filt,features = "LYZ",pt.size = 0)



## ----clust_resoFix1,include=TRUE,eval=T, fig.height=4,fig.width=7, warning=FALSE, message=FALSE----

seu_filt <- FindClusters(seu_filt, resolution = 0.6)


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Optimized Clusters")


## ----include=TRUE,eval=T, fig.height=4,fig.width=7----------------------------
DimPlot(seu_filt, group.by = "cluster_byDefault",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Default Clusters")


## ----clust_resoFix2,include=TRUE,eval=T, fig.height=4,fig.width=7, warning=FALSE, message=FALSE----
seu_filt <- FindClusters(seu_filt, resolution = 0.4)


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Optimized Clusters")


## ----include=TRUE,eval=T, fig.height=4,fig.width=7----------------------------
DimPlot(seu_filt, group.by = "cluster_byDefault",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Default Clusters")


## ----markGene_cal,include=TRUE,eval=F-----------------------------------------
# markers <- FindAllMarkers(seu_filt, only.pos = TRUE,
#                           min.pct = 0.25, logfc.threshold = 0.25)


## ----eval=F, echo=F-----------------------------------------------------------
# save(markers, file="data/pbmc8k_markers_updated.RData")


## ----echo=F, eval=T-----------------------------------------------------------
load("data/pbmc8k_markers_updated.RData")


## -----------------------------------------------------------------------------
head(markers)


## ----top2Mark,include=TRUE,eval=T---------------------------------------------
top_genes <- markers %>% group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
head(top_genes)



## ----fig.height=4,fig.width=7-------------------------------------------------

DoHeatmap(seu_filt, features = top_genes$gene) + NoLegend()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Evaluate known marker genes expression

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html>

---
"
  )
}else{
  cat("# Evaluate known marker genes expression


---
"
  )

}



## ----resEval_knownMarker_allMark,include=TRUE,eval=T,fig.height=4,fig.width=7----
known_marker <- c("IL7R","CCR7")
FeaturePlot(seu_filt, features = known_marker)


## ----include=TRUE,eval=T,fig.height=4,fig.width=7-----------------------------
known_marker <- c("GNLY","NKG7")
FeaturePlot(seu_filt, features = known_marker)


## ----resEval_knowMarker_heatmap,include=TRUE,eval=T---------------------------
known_marker <- c("IL7R","CCR7","S100A4","CD8A",
                                   "MS4A1",
                                   "CD14","LYZ","FCGR3A","MS4A7",
                                   "FCER1A","CST3",
                                   "GNLY","NKG7","PPBP")

mat <- GetAssayData(seu_filt,assay = "RNA",slot = "data")
mat <- mat[known_marker,]
mat <- as.matrix(mat)



## -----------------------------------------------------------------------------
clust <- unique(seu_filt$seurat_clusters)
clust <- as.character(clust)

avgExp_byClust <- lapply(clust,function(clust, seu, known){
  sub <- subset(seu, subset=seurat_clusters==clust)
  mat <- GetAssayData(sub,assay="RNA",slot = "data")
  mat <- mat[known,]
  mat <- as.matrix(mat)
  avg <- rowMeans(mat)
  return(avg)}, seu_filt, known_marker)

names(avgExp_byClust) <- paste0("C",clust)



## ----fig.height=4,fig.width=7-------------------------------------------------
library(pheatmap)
avgExp_mat <- do.call(cbind, avgExp_byClust)
pheatmap(avgExp_mat,scale = "row")


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(avgExp_mat, avgExp_byClust, clust, mat, marker)
gc()


## ----sec2_resEval_cellTypeAsign,include=TRUE,eval=T---------------------------

seu_filt[["cellType_byClust"]] <- NA
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(0,5)] <- "Naive CD4+ T cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(3)] <- "Memory CD4+ T cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(4)] <- "CD8+ T cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(6)] <- "NK cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(7,10)] <- "DC"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(8)] <- "FCGR3A+ Monocytes"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(2)] <- "B cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(1)] <- "CD14+ Monocytes"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(11)] <- "Platelets"
seu_filt$cellType_byClust[is.na(seu_filt$cellType_byClust)] <- "Unspecified"


## -----------------------------------------------------------------------------

table(seu_filt$cellType_byClust)


## ----fig.height=4,fig.width=7-------------------------------------------------

DimPlot(seu_filt, group.by = c("seurat_clusters","cellType_byClust"),label = F,pt.size = 0.2)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Loupe Browser

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("#  Loupe Browser

---
"    
  )
  
}



## ----eval=F-------------------------------------------------------------------
# remotes::install_github("10xGenomics/loupeR")
# loupeR::setup()


## ----eval=F-------------------------------------------------------------------
# library(loupeR)
# create_loupe_from_seurat(seu_filt,
#                          output_dir = "loupe",
#                          output_name = "seu_filt")
# 

