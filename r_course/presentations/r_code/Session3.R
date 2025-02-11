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
knitr::opts_chunk$set(echo = TRUE, tidy = T, fig.height=4, fig.width=7, warning = F, message=F)


## ----sec3_loadPack,include=F,eval=FALSE---------------------------------------
# library(Seurat)
# library(scran)
# library(scater)
# library(ggplot2)
# library(pheatmap)
# library(TSCAN)
# library(SoupX)
# library(DropletUtils)
# library(scuttle)


## ----sec_mergeData_funcUsed_dataProc,include=TRUE-----------------------------
data_proc <- function(seu){
  seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, select.method="vst", nfeatures=2000)
  return(seu)}


## ----sec3_mergeData_funcUsed_quickClust,include=TRUE--------------------------
quick_clust <- function(seu){
  set.seed(42)
  seu <- ScaleData(seu, verbose=FALSE)
  seu <- RunPCA(seu, npcs=30, verbose=FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, verbose=FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:10, verbose=FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose=FALSE)
  return(seu)}


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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



## -----------------------------------------------------------------------------
my_seu_list <- readRDS("data/to_integrate.rds")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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

seu_merge <- merge(my_seu_list[[1]], my_seu_list[2:4],
                            add.cell.ids = c("AD2b", "AD4", "C1", "C3"), project = "Merge")
head(seu_merge,4)



## ----sec3_mergeData_woCorr_cluster,include=TRUE-------------------------------
seu_merge <- data_proc(seu_merge)
seu_merge <- quick_clust(seu_merge)


## -----------------------------------------------------------------------------
DimPlot(seu_merge, group.by = "sample_id", pt.size = 0.2)





## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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



feats <- SelectIntegrationFeatures(my_seu_list)

my_seu_list_rpca <- lapply(my_seu_list, function(seu, feats){
  seu <- ScaleData(seu, features=feats, verbose=FALSE)
  seu <- RunPCA(seu, features=feats, verbose=FALSE)
  return(seu)}, feats)


## ----sec3_mergeData_RPCA_int,include=TRUE,eval=T------------------------------
anchors <- FindIntegrationAnchors(my_seu_list_rpca, anchor.features = feats, reduction = "rpca")

my_seu_merge_rpca <- IntegrateData(anchorset = anchors)

my_seu_merge_rpca


## ----sec3_mergeData_RPCA_cluster,include=TRUE,eval=T--------------------------

my_seu_merge_rpca <- ScaleData(my_seu_merge_rpca)
my_seu_merge_rpca <- quick_clust(my_seu_merge_rpca)


## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "sample_id", pt.size = 0.2)



## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "seurat_clusters", pt.size = 0.2)



## -----------------------------------------------------------------------------
library(pheatmap)
tbl <- table(my_seu_merge_rpca$sample_id, my_seu_merge_rpca$seurat_clusters)
pheatmap(tbl, scale="row")


## -----------------------------------------------------------------------------
FeaturePlot(my_seu_merge_rpca, "MOBP", pt.size = 0.2)



## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "paper_annot", pt.size = 0.2)



## -----------------------------------------------------------------------------
library(pheatmap)

tbl <- table(my_seu_merge_rpca$seurat_clusters, my_seu_merge_rpca$paper_annot)
pheatmap(tbl, scale = "row")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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

seu_merge <- merge(my_seu_list[[1]], my_seu_list[2:4],
                            add.cell.ids = c( "C1", "C3","AD2b", "AD4"), project = "Merge")
seu_merge <- data_proc(seu_merge)
seu_merge <- ScaleData(seu_merge)
seu_merge <- RunPCA(seu_merge)
seu_merge <- RunUMAP(seu_merge, reduction = "pca", dims = 1:10, reduction.name = "umap")


## -----------------------------------------------------------------------------
DimPlot(seu_merge, group.by = "sample_id")


## ----sec3_mergedata_Harmony_merge,include=TRUE,eval=T-------------------------
library(harmony)
seu_merge_harmony <- RunHarmony(seu_merge, group.by.vars= "sample_id", assay.use= "RNA")
seu_merge_harmony <- RunUMAP(seu_merge_harmony, reduction = "harmony", dims = 1:10, reduction.name = "umap_harmony")
seu_merge_harmony <- FindNeighbors(seu_merge_harmony, reduction = "harmony")
seu_merge_harmony <- FindClusters(seu_merge_harmony) 
    



## -----------------------------------------------------------------------------
DimPlot(seu_merge_harmony, group.by = "sample_id", reduction = "umap_harmony", pt.size = 0.2)



## -----------------------------------------------------------------------------
DimPlot(seu_merge_harmony, group.by = "seurat_clusters", reduction = "umap_harmony",  pt.size = 0.2)



## -----------------------------------------------------------------------------
library(pheatmap)
tbl <- table(seu_merge_harmony$sample_id, seu_merge_harmony$seurat_clusters)
pheatmap(tbl, scale="row")


## -----------------------------------------------------------------------------
FeaturePlot(seu_merge_harmony, "MOBP", pt.size = 0.2)



## -----------------------------------------------------------------------------
DimPlot(seu_merge_harmony, group.by = "paper_annot", pt.size = 0.2)



## -----------------------------------------------------------------------------

tbl <- table(seu_merge_harmony$seurat_clusters, seu_merge_harmony$paper_annot)
pheatmap(tbl, scale = "row")



## ----echo=F, eval=F-----------------------------------------------------------
# 
# saveRDS(my_seu_merge_rpca,"~/Desktop/integrated.rds" )
# saveRDS(my_seu_merge_rpca, "scRNASeq/inst/extdata/data/integrated.rds")
# 


## ----echo=F-------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "sample_id") + ggtitle("rPCA")



## ----echo=F-------------------------------------------------------------------
DimPlot(seu_merge_harmony, group.by = "sample_id", reduction="umap_harmony") + ggtitle("harmony")


## ----echo=F-------------------------------------------------------------------
tbl <- table(my_seu_merge_rpca$seurat_clusters, my_seu_merge_rpca$paper_annot)
pheatmap(tbl, scale = "row", main="rPCA", treeheight_row = 0, , treeheight_col = 0)


## ----echo=F-------------------------------------------------------------------
tbl <- table(seu_merge_harmony$seurat_clusters, seu_merge_harmony$paper_annot)
pheatmap(tbl, main="Harmony", treeheight_row = 0, , treeheight_col = 0)


## ----echo=F-------------------------------------------------------------------
  tbl <- table(my_seu_merge_rpca$sample_id, my_seu_merge_rpca$seurat_clusters)
pheatmap(tbl, scale = "row", main="rPCA", treeheight_row = 0, , treeheight_col = 0)



## ----echo=F-------------------------------------------------------------------
tbl <- table(seu_merge_harmony$sample_id, seu_merge_harmony$seurat_clusters)
pheatmap(tbl, scale = "row", main="Harmony", treeheight_row = 0, , treeheight_col = 0)


## ----echo=F, warning=FALSE, message=F-----------------------------------------

rm(seu_merge_harmony, my_seu_list_rpca,my_seu_list,seu_merge)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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



## -----------------------------------------------------------------------------

DimPlot(my_seu_merge_rpca, group.by = "paper_annot")


## -----------------------------------------------------------------------------

FeaturePlot(my_seu_merge_rpca, c("MOBP","MAG"))


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Annotation with SingleR

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Annotation with SingleR

---
"    
  )
  
}



## -----------------------------------------------------------------------------
library(celldex)
hpcad <- HumanPrimaryCellAtlasData()



## -----------------------------------------------------------------------------
hpcad 



## -----------------------------------------------------------------------------
my_seu_merge_rpca



## -----------------------------------------------------------------------------
my_seu_merge_rpca_mat <- GetAssayData(my_seu_merge_rpca)

my_seu_merge_rpca_mat[1:5,1:5]


## -----------------------------------------------------------------------------
colData(hpcad)



## -----------------------------------------------------------------------------
library(SingleR)

pred_res <- SingleR(ref = hpcad, test = my_seu_merge_rpca_mat, labels = hpcad$label.main)



## -----------------------------------------------------------------------------
head(pred_res,2)


## -----------------------------------------------------------------------------
mat <- as.matrix(pred_res$scores)
rownames(mat) <- rownames(pred_res)
pheatmap::pheatmap(mat, scale = "row", show_rownames = FALSE, fontsize_col = 5)


## -----------------------------------------------------------------------------
my_seu_merge_rpca$hpcad_singleR_labels <- pred_res$labels

summ_table <- table(my_seu_merge_rpca$hpcad_singleR_labels, my_seu_merge_rpca$paper_annot)

pheatmap(summ_table, scale="column", fontsize_row = 5)


## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "hpcad_singleR_labels")



## ----eval=F-------------------------------------------------------------------
# abm <- readRDS("data/abm.rds")
# 


## ----eval=F,echo=TRUE---------------------------------------------------------
# head(abm)
# 


## ----echo=F-------------------------------------------------------------------
out<-read.csv("data/out.txt", row.names=1)
out


## ----eval=F, echo=T-----------------------------------------------------------
# 
# pred_res2 <- SingleR(ref = GetAssayData(abm), test = my_seu_merge_rpca_mat, labels = abm$class_label)


## ----eval=F, echo=F-----------------------------------------------------------
# saveRDS(pred_res2,file = "data/annotate_df1.rds")


## ----eval=T, echo=F-----------------------------------------------------------
pred_res2 <- readRDS(file = "data/annotate_df1.rds")
my_seu_merge_rpca <- readRDS(file = "data/annotated.rds")


## -----------------------------------------------------------------------------
head(pred_res2,2)


## -----------------------------------------------------------------------------
mat <- as.matrix(pred_res2$scores)
rownames(mat) <- rownames(pred_res2)
pheatmap::pheatmap(mat, scale = "row", show_rownames = FALSE, fontsize_col = 5)


## ----eval=F-------------------------------------------------------------------
# 
# my_seu_merge_rpca$abm_singleR_labels <- pred_res2$labels


## -----------------------------------------------------------------------------
summ_table <- table(my_seu_merge_rpca$abm_singleR_labels, my_seu_merge_rpca$paper_annot)

pheatmap(summ_table, scale="column", fontsize_row = 5)


## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "abm_singleR_labels")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Annotation with Seurat

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Annotation with Seurat

---
"    
  )
  
}



## ----sec3_ctAnno_Seurat_gatherData,include=TRUE,eval=F------------------------
# 
# transfer_list <- list("ref"=abm,"query"=my_seu_merge_rpca)
# feats <- SelectIntegrationFeatures(transfer_list)
# 
# transfer_list <- lapply(transfer_list,function(seu,feats){
#   seu <- ScaleData(seu,features=feats,verbose=FALSE)
#   seu <- RunPCA(seu,features=feats,verbose=FALSE)
#   return(seu)}, feats)
# 


## ----sec3_ctAnno_Seurat_tranferAnno,include=TRUE,eval=F-----------------------
# anchors <- FindTransferAnchors(reference = transfer_list$ref,
#                                query=transfer_list$query,
#                                dims=1:30, reference.reduction="pca")
# pred_res3 <- TransferData(anchorset = anchors, refdata=transfer_list$ref$class_label)


## ----eval=F, echo=F-----------------------------------------------------------
# saveRDS(pred_res3,file = "data/annotate_df2.rds")


## ----eval=T, echo=F-----------------------------------------------------------
pred_res3 <- readRDS(file = "data/annotate_df2.rds")



## -----------------------------------------------------------------------------
head(pred_res3,2)


## -----------------------------------------------------------------------------
mat <- as.matrix(pred_res3[,-c(1,5)])
colnames(mat) <- gsub("prediction.score.","",colnames(mat))
pheatmap(mat,scale = "row",show_rownames = FALSE)


## ----eval=F-------------------------------------------------------------------
# 
# my_seu_merge_rpca$abm_seurat_labels <- pred_res3$predicted.id
# 


## -----------------------------------------------------------------------------

summ_table <- table(my_seu_merge_rpca$abm_seurat_labels, my_seu_merge_rpca$paper_annot)

pheatmap(summ_table, scale="column", fontsize_row = 5)


## -----------------------------------------------------------------------------
DimPlot(my_seu_merge_rpca, group.by = "abm_seurat_labels")



## -----------------------------------------------------------------------------

table(my_seu_merge_rpca$abm_seurat_labels == my_seu_merge_rpca$abm_singleR_labels)
tbl <- table(my_seu_merge_rpca$abm_seurat_labels,my_seu_merge_rpca$abm_singleR_labels)
pheatmap::pheatmap(tbl,scale = "row")

