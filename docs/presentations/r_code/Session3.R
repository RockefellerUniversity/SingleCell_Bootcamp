params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(DropletUtils))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(scuttle))
suppressPackageStartupMessages(require(scDblFinder))
knitr::opts_chunk$set(echo = TRUE, tidy = T, fig.height=4, fig.width=7)


## ----sec3_loadPack,include=TRUE,eval=FALSE------------------------------------
## library(Seurat)
## library(scran)
## library(scater)
## library(SeuratData)
## library(batchelor)
## library(ggplot2)
## library(pheatmap)
## library(slingshot)
## library(TSCAN)
## library(SoupX)
## library(DropletUtils)
## library(scuttle)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Merging Datasets

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Merging Datasets

---
"    
  )
  
}



## ---- eval=F------------------------------------------------------------------
## devtools::install_github('satijalab/seurat-data')


## -----------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
InstallData("ifnb")
LoadData("ifnb")
head(ifnb,2)


## ----sec3_mergeData_fetchEG,include=TRUE, eval=F------------------------------
## table(ifnb$stim)
## ifnb_list <- Seurat::SplitObject(ifnb, split.by="stim")
## ifnb_list


## ---- echo=F, eval=F----------------------------------------------------------
## save("ifnb_list",
##      file = "data/seuOBJ_IFNB_splitByStim.RData")


## -----------------------------------------------------------------------------
load("data/seuOBJ_IFNB_splitByStim.RData")



## ----sec_mergeData_funcUsed_dataProc,include=TRUE-----------------------------
data_proc <- function(seu){
  seu <- NormalizeData(seu,normalization.method="LogNormalize",scale.factor=10000)
  seu <- FindVariableFeatures(seu,select.method="vst",nfeatures=2000)
  return(seu)}


## ----sec3_mergeData_funcUsed_quickClust,include=TRUE--------------------------
quick_clust <- function(seu){
  set.seed(1001)
  seu <- ScaleData(seu,verbose=FALSE)
  seu <- RunPCA(seu,npcs=30,verbose=FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:10,verbose=FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:10,verbose=FALSE)
  seu <- FindClusters(seu, resolution = 0.5,verbose=FALSE)
  return(seu)}


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Simple Merge

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Simple Merge

---
"    
  )
  
}



## ----mergeData_woCor_merged,include=TRUE,eval=T-------------------------------

ifnb_merge <- merge(ifnb_list$CTRL, ifnb_list$STIM,
                            add.cell.ids = c("CTRL","STIM"), project = "ifnb_seuMerge")
head(ifnb_merge,2)


## ----sec3_mergeData_woCorr_cluster,include=TRUE-------------------------------
ifnb_merge <- data_proc(ifnb_merge)
ifnb_merge <- quick_clust(ifnb_merge)


## -----------------------------------------------------------------------------
DimPlot(ifnb_merge,group.by = "stim", pt.size = 0.2)


## ----sec3_mergeData_woCorr_eval,include=TRUE----------------------------------

DimPlot(ifnb_merge, group.by = "seurat_annotations", pt.size = 0.2,split.by = "stim")



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Merge with reciprocal PCA

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Merge with reciprocal PCA

---
"    
  )
  
}



## ----sec3_mergeData_RPCA_prep,include=TRUE,eval=T-----------------------------

ifnb_list_rpca <- lapply(ifnb_list, data_proc)

feats <- SelectIntegrationFeatures(ifnb_list_rpca)

ifnb_list <- lapply(ifnb_list,function(seu,feats){
  seu <- ScaleData(seu,features=feats,verbose=FALSE)
  seu <- RunPCA(seu,features=feats,verbose=FALSE)
  return(seu)},feats)


## ----sec3_mergeData_RPCA_int,include=TRUE,eval=T------------------------------
anchors <- FindIntegrationAnchors(ifnb_list, anchor.features = feats, reduction = "rpca")

ifnb_merge <- IntegrateData(anchorset = anchors)

ifnb_merge


## ----sec3_mergeData_RPCA_cluster,include=TRUE,eval=F--------------------------
## 
## ifnb_merge <- ScaleData(ifnb_merge)
## 
## ifnb_merge <- quick_clust(ifnb_merge)
## 
## DimPlot(ifnb_merge,group.by = "stim",pt.size = 0.2)


## ----sec3_mergeData_RPCA_cluster2,include=TRUE,eval=T, echo=F-----------------

ifnb_merge <- ScaleData(ifnb_merge)

ifnb_merge <- quick_clust(ifnb_merge)

DimPlot(ifnb_merge,group.by = "stim",pt.size = 0.2)


## -----------------------------------------------------------------------------
tbl <- table(ifnb_merge$stim, ifnb_merge$seurat_clusters)

barplot(tbl,beside = T, main= "Cell numbers in each cluster of each group")


## ----sec3_mergeData_RPCA_eval,include=TRUE,eval=T-----------------------------

DimPlot(ifnb_merge, group.by = "seurat_annotations", split.by= "stim", pt.size = 0.2)


## -----------------------------------------------------------------------------
library(pheatmap)

tbl <- table(ifnb_merge$seurat_clusters,ifnb_merge$seurat_annotations)
pheatmap(tbl, scale = "column")



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Merge data with MNN correction

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Merge data with MNN correction

---
"    
  )
  
}



## ----sec3_mergeData_MNN_prep,include=TRUE,eval=T------------------------------

sce_list <- lapply(ifnb_list,function(seu){
  sce <- as.SingleCellExperiment(seu,assay="RNA")
  rowData(sce)$SYMBOL <- rownames(sce)
  return(sce)})


sce_list 


## -----------------------------------------------------------------------------
library(scran)
dec_list <- lapply(sce_list, modelGeneVar)

hvgc_list <-lapply(sce_list, getTopHVGs ,prop=0.1)


## -----------------------------------------------------------------------------
universe <- intersect(rownames(sce_list$CTRL),rownames(sce_list$STIM))
sce_list <- lapply(sce_list,function(sce,universe){
  sce <- sce[universe,];return(sce)},universe)
dec_list <- lapply(dec_list,function(dec,universe){
  dec <- dec[universe,];return(dec)},universe)

combined_dec <- combineVar(dec_list$CTRL, dec_list$STIM)
chosen_hvgs <- combined_dec$bio > 0


## ----sec3_mergeData_MNN_cor,include=TRUE,eval=T-------------------------------
library(batchelor)
mnn_res <- fastMNN("CTRL"=sce_list$CTRL, "STIM"=sce_list$STIM,
                   d=50, 
                   k=20,
                   subset.row=chosen_hvgs)
mnn_res


## -----------------------------------------------------------------------------
table(mnn_res$batch)


## -----------------------------------------------------------------------------
reducedDim(mnn_res,"corrected")[1:2,]


## -----------------------------------------------------------------------------
assay(mnn_res,"reconstructed")[1:2,1:5]


## ----sec3_dataMerge_MNN_cluster,include=TRUE,eval=T---------------------------
library(scater)
set.seed(1001)
mnn_res <- runUMAP(mnn_res, dimred="corrected")
mnn_res$batch <- factor(mnn_res$batch)
plotUMAP(mnn_res, colour_by="batch")


## -----------------------------------------------------------------------------
snn.gr <- buildSNNGraph(mnn_res,use.dimred="corrected")
cluster_mnn <- igraph::cluster_louvain(snn.gr)$membership
table(cluster_mnn)


## -----------------------------------------------------------------------------
mnn_res$cluster <- factor(cluster_mnn)
plotUMAP(mnn_res,colour_by="cluster")


## -----------------------------------------------------------------------------

tbl <- table(mnn_res$batch,mnn_res$cluster)
barplot(tbl,beside = T, main= "Cell numbers in each cluster of each group")


## ----sec3_dataMerge_MNN_eval,include=TRUE,eval=T------------------------------

cellType <- lapply(sce_list,function(x){
  res <- setNames(as.character(colData(x)$seurat_annotations),colnames(x))
  return(res)})
cell_type <- c(cellType$CTRL,cellType$STIM)
mnn_res$cell_type <- cell_type[match(rownames(colData(mnn_res)),names(cell_type))]
mnn_res$cell_type <- factor(mnn_res$cell_type)


## -----------------------------------------------------------------------------
plotUMAP(mnn_res,colour_by="cell_type")


## -----------------------------------------------------------------------------

plotUMAP(mnn_res,colour_by="cell_type")


## -----------------------------------------------------------------------------
tbl <- table(mnn_res$cluster, mnn_res$cell_type)
pheatmap(tbl,scale = "column")



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Merge data with Harmony

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Merge data with Harmony

---
"    
  )
  
}



## ----sec3_mergedata_Harmony_samplePrep,include=TRUE,eval=T--------------------

seu_obj <- merge(ifnb_list$CTRL, ifnb_list$STIM)
seu_obj <- data_proc(seu_obj)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj)
seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims = 1:10, reduction.name = "umap")


## -----------------------------------------------------------------------------
DimPlot(seu_obj)


## ----sec3_mergedata_Harmony_merge,include=TRUE,eval=FALSE---------------------
## library(harmony)
## seu_obj <- RunHarmony(seu_obj, group.by.vars = "stim", assay.use= "RNA")
## seu_obj <- RunUMAP(seu_obj, reduction = "harmony", dims = 1:10, reduction.name = "umap_harmony")
## DimPlot(seu_obj, reduction = "umap_harmony")


## ---- echo=F, eval=T, warning=F, message=FALSE, include=F---------------------

rm(seu_obj, ifnb, ifnb_list, ifnb_merge, ifnb_list_rpca, feats, anchors, sce_list, dec_list, hvgc_list, combined_dec, chosen_hvgs, mnn_res, snn.gr, cluster_mnn ,cellType )
gc()


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell type annotation

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell type annotation

---
"    
  )
  
}



## ----sec3_ctAnno_prepDemo_loadData,include=TRUE,eval=T------------------------
library(Seurat)
library(SeuratData)
InstallData("panc8")
data("panc8")
head(panc8,2)


## -----------------------------------------------------------------------------
seu_list <- SplitObject(panc8, split.by = "tech")
names(seu_list) <- c("celseq1", "celseq2", "smartseq", "fluidigmc1", "indrop")


## ----sec3_ctAnno_prepDemo_buildREF,include=TRUE,eval=T------------------------

seu_list <- lapply(seu_list, data_proc)
ref_list <- seu_list[c("celseq1","celseq2")]
feats <- SelectIntegrationFeatures(ref_list)
ref_list <- lapply(ref_list,function(seu,feats){
  seu <- ScaleData(seu,features=feats,verbose=FALSE)
  seu <- RunPCA(seu,features=feats,verbose=FALSE)
  return(seu)},feats)


## ----sec3_ctAnno_Seurat_RPCA,include=TRUE,eval=T------------------------------

anchors <- FindIntegrationAnchors(ref_list,anchor.features = feats,reduction = "rpca")
ref_panc <- IntegrateData(anchorset = anchors)
ref_panc <- ScaleData(ref_panc)
ref_panc <- quick_clust(ref_panc)


## -----------------------------------------------------------------------------
DimPlot(ref_panc,group.by = "seurat_clusters",split.by = "tech",pt.size = 0.2,label=TRUE)


## ----sec3_ctAnno_Seurat_gatherData,include=TRUE,eval=T------------------------

panc_list <- list("ref"=ref_panc,"query"=seu_list$smartseq)
feats <- SelectIntegrationFeatures(panc_list)
panc_list <- lapply(panc_list,function(seu,feats){
  seu <- ScaleData(seu,features=feats,verbose=FALSE)
  seu <- RunPCA(seu,features=feats,verbose=FALSE)
  return(seu)},feats)


## ----sec3_ctAnno_Seurat_tranferAnno,include=TRUE,eval=T-----------------------
anchors <- FindTransferAnchors(reference = panc_list$ref,
                               query=panc_list$query,
                               dims=1:30,reference.reduction="pca")
pred_res <- TransferData(anchorset = anchors,refdata=panc_list$ref$celltype)
head(pred_res,2)


## -----------------------------------------------------------------------------
mat <- as.matrix(pred_res[,-c(1,15)])
colnames(mat) <- gsub("prediction.score.","",colnames(mat))
pheatmap::pheatmap(mat,scale = "row",show_rownames = FALSE)


## ----sec3_ctAnno_Seurat_pres,include=TRUE,eval=T------------------------------

pred_cellType <- setNames(pred_res$predicted.id, rownames(pred_res))
panc_list$query[["cellType_predBySeurat"]] <- pred_cellType[match(Cells(panc_list$query),
                                                                  names(pred_cellType))]

head(panc_list$query,2)

table(panc_list$query$cellType_predBySeurat)


## ----sec3_ctAnno_SingleR_pred,include=TRUE,eval=T-----------------------------
sce_list <- lapply(panc_list, function(seu){
  sce <- as.SingleCellExperiment(seu, assay="RNA")
  return(sce)})



## -----------------------------------------------------------------------------
library(SingleR)
pred_res <- SingleR(ref = sce_list$ref, test = sce_list$query, labels = sce_list$ref$celltype)
head(pred_res,2)


## -----------------------------------------------------------------------------
mat <- as.matrix(pred_res$scores)
rownames(mat) <- rownames(pred_res)
pheatmap::pheatmap(mat, scale = "row", show_rownames = FALSE)


## ----sec3_ctAnno_SingleR_import,include=TRUE,eval=T---------------------------
cell_type <- setNames(pred_res$pruned.labels,rownames(pred_res))
panc_list$query$cellType_predBySingleR <- cell_type[match(Cells(panc_list$query),names(cell_type))]
head(panc_list$query,2)



## -----------------------------------------------------------------------------
table(panc_list$query$cellType_predBySeurat == panc_list$query$celltype)


## -----------------------------------------------------------------------------
table(panc_list$query$cellType_predBySingleR == panc_list$query$celltype)


## -----------------------------------------------------------------------------

table(panc_list$query$cellType_predBySeurat == panc_list$query$cellType_predBySingleR)
tbl <- table(panc_list$query$cellType_predBySeurat,panc_list$query$cellType_predBySingleR)
pheatmap::pheatmap(tbl,scale = "row")


## ---- echo=F, eval=T, warning=F, message=FALSE, include=F---------------------

rm(panc8, seu_list, ref_list, feats, anchors, ref_panc, panc_list, seu, pred_res, mat, pred_cellType, sce_list, cell_type)
gc()


## ---- echo=F, eval=T, warning=F, message=FALSE, include=F---------------------

rm(sce.sling2)
rm(res, embedded, embedded_curve,sce2, pseudo.paths, avg_pseudoTime, curve, sce.sling, sce, tbl, panc_list, cell_type, mat, pred_res,sce_list,pred_cellType,anchors,ref_panc,feats,ref_list,seu_list,panc8, seu_obj,mnn_res,cell_type,snn.gr,cluster_mnn,universe,sce_list,sce,dec_list,dec,combined_dec,chosen_hvgs,hvgc_list,anchors,ifnb_merge,ifnb_list,feats,ifnb_list_rpca ,ifnb)
gc()



