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



## ----eval=F-------------------------------------------------------------------
## remotes::install_github('satijalab/seurat-data')
## 


## ----eval=F-------------------------------------------------------------------
## library(Seurat)
## library(SeuratData)
## InstallData("ifnb")
## LoadData("ifnb")
## head(ifnb,2)


## ----echo=F, eval=F-----------------------------------------------------------
## write.csv(head(ifnb,2),"data/ifnb_head.csv")
## 


## ----echo=F, eval=T-----------------------------------------------------------
library(Seurat)
ifnb_test <- read.csv("data/ifnb_head.csv",row.names=1)
ifnb_test


## ----echo=F, eval=F-----------------------------------------------------------
## 
## # Download: https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz
## install.packages("~/Downloads/ifnb.SeuratData_3.0.2.tar.gz", repos = NULL, type = "source")
## library(ifnb.SeuratData)
## ifnb<-LoadData("ifnb")
## 


## ----sec3_mergeData_fetchEG,include=TRUE, eval=F------------------------------
## table(ifnb$stim)
## ifnb_list <- Seurat::SplitObject(ifnb, split.by="stim")
## ifnb_list


## ----echo=F, eval=F-----------------------------------------------------------
## save("ifnb_list",
##      file = "data/seuOBJ_IFNB_splitByStim.RData")


## ----eval=F-------------------------------------------------------------------
## 
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
  set.seed(42)
  seu <- ScaleData(seu,verbose=FALSE)
  seu <- RunPCA(seu,npcs=30,verbose=FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:10,verbose=FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:10,verbose=FALSE)
  seu <- FindClusters(seu, resolution = 0.5,verbose=FALSE)
  return(seu)}


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

ifnb_merge <- merge(ifnb_list$CTRL, ifnb_list$STIM,
                            add.cell.ids = c("CTRL","STIM"), project = "ifnb_seuMerge")
head(ifnb_merge,2)


## ----eval=F, echo=F-----------------------------------------------------------
## remotes::install_version("Matrix", version = "1.6-1.1")


## ----sec3_mergeData_woCorr_cluster,include=TRUE-------------------------------
ifnb_merge <- data_proc(ifnb_merge)
ifnb_merge <- quick_clust(ifnb_merge)


## -----------------------------------------------------------------------------
DimPlot(ifnb_merge,group.by = "stim", pt.size = 0.2)


## ----sec3_mergeData_woCorr_eval,include=TRUE----------------------------------

DimPlot(ifnb_merge, group.by = "seurat_annotations", pt.size = 0.2,split.by = "stim")



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



## ----sec3_mergeData_RPCA_prep,include=TRUE,eval=F-----------------------------
## 
## ifnb_list_rpca <- lapply(ifnb_list, data_proc)
## 
## feats <- SelectIntegrationFeatures(ifnb_list_rpca)
## 
## ifnb_list <- lapply(ifnb_list,function(seu,feats){
##   seu <- ScaleData(seu,features=feats,verbose=FALSE)
##   seu <- RunPCA(seu,features=feats,verbose=FALSE)
##   return(seu)},feats)


## ----eval=F, echo=F-----------------------------------------------------------
## saveRDS(ifnb_list, "data/ifnb_list.rds")
## saveRDS(feats, "data/feats.rds")
## 


## ----eval=T, echo=F-----------------------------------------------------------
ifnb_list <- readRDS("data/ifnb_list.rds")
feats <- readRDS("data/feats.rds")



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



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
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
tbl <- table(mnn_res$cluster, mnn_res$cell_type)
pheatmap(tbl,scale = "column")



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


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu_obj, ifnb, ifnb_list, ifnb_merge, ifnb_list_rpca, anchors, sce_list, dec_list, hvgc_list, combined_dec, chosen_hvgs, mnn_res, snn.gr, cluster_mnn ,cellType )
gc()


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



## ----sec3_ctAnno_prepDemo_loadData,include=TRUE,eval=F------------------------
## library(Seurat)
## library(SeuratData)
## InstallData("panc8")
## data("panc8")
## 


## ----echo=F, eval=F-----------------------------------------------------------
## 
## # Download: https://seurat.nygenome.org/src/contrib/panc8.SeuratData_3.0.2.tar.gz
## install.packages("~/Downloads/panc8.SeuratData_3.0.2.tar.gz", repos = NULL, type = "source")
## library(panc8.SeuratData)
## panc8 <- LoadData("panc8")
## UpdateSeuratObject(panc8)


## ----eval=F-------------------------------------------------------------------
## panc8 <- UpdateSeuratObject(panc8)
## seu_list <- SplitObject(panc8, split.by = "tech")
## names(seu_list) <- c("celseq1", "celseq2", "smartseq", "fluidigmc1", "indrop")


## ----eval=F, echo=F-----------------------------------------------------------
## save(seu_list, file="data/panc8.RData")


## -----------------------------------------------------------------------------
load("data/panc8.RData")


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


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(panc8, seu_list, ref_list, anchors, ref_panc, panc_list, seu, pred_res, mat, pred_cellType, sce_list, cell_type)
gc()


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



## ----eval=T-------------------------------------------------------------------
download.file("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz","10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz")


## ----eval=T-------------------------------------------------------------------
untar("10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz")



## ----sec3_dropProc_kneePlot_prep,include=TRUE, eval=T-------------------------
library(DropletUtils)
library(Seurat)

raw_mtx <- Seurat::Read10X("raw_feature_bc_matrix")
dim(raw_mtx)



## -----------------------------------------------------------------------------
bcrank <- barcodeRanks(raw_mtx)
head(bcrank,2)


## ----sec3_dropProc_kennPlot_press,include=TRUE,eval=T-------------------------

uniq <- !duplicated(bcrank$rank)
bcrank <- bcrank[uniq,]

knee <- metadata(bcrank)$knee
message(paste0("knee point: ",knee))
inflection <- metadata(bcrank)$inflection
message(paste0("inflection point: ",inflection))


## ----eval=F-------------------------------------------------------------------
## ggplot(as.data.frame(bcrank),aes(x=rank,y=total)) + geom_point()+
##   geom_hline(yintercept = knee, linetype="dashed",color="blue")+
##   geom_hline(yintercept = inflection, linetype="dashed",color="green")+
##   scale_x_continuous(trans = "log10")+
##   scale_y_continuous(trans = "log10")+
##   labs(x="cell barcode ranked by counts",
##        y="UMI counts of each cell barcode")+
##   theme_classic()


## ----eval=T, echo=F-----------------------------------------------------------
ggplot(as.data.frame(bcrank),aes(x=rank,y=total)) + geom_point()+
  geom_hline(yintercept = knee, linetype="dashed",color="blue")+
  geom_hline(yintercept = inflection, linetype="dashed",color="green")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  labs(x="cell barcode ranked by counts",
       y="UMI counts of each cell barcode")+
  theme_classic()


## ----sec3_dropProc_remEmtpy_cal, eval=TRUE------------------------------------

e.out <- emptyDrops(raw_mtx)
e.out <- e.out[order(e.out$FDR),]
head(e.out)



## ----eval=F-------------------------------------------------------------------
## 
## ggplot(as.data.frame(e.out),aes(x=Total))+geom_histogram()+
##   geom_vline(xintercept = knee, color="red",linetype="dashed")+
##   labs(x="UMI counts per cell",y="Frequency")+
##   scale_y_continuous(trans = "log10")+
##   scale_x_continuous(trans = "log10")+
##   theme_classic()


## ----eval=T, echo=F-----------------------------------------------------------

ggplot(as.data.frame(e.out),aes(x=Total))+geom_histogram()+
  geom_vline(xintercept = knee, color="red",linetype="dashed")+
  labs(x="UMI counts per cell",y="Frequency")+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  theme_classic()


## -----------------------------------------------------------------------------
table(e.out$FDR < 0.001)
cell_sel <- rownames(e.out)[e.out$FDR < 0.001]
filt_mtx <- raw_mtx[,colnames(raw_mtx) %in% cell_sel]


## ----sec3_dropProc_ambIDENT_DorpUtils,include=TRUE,eval=T---------------------
amb <- metadata(e.out)$ambient[,1]
head(amb)


## ----eval=F, echo=F-----------------------------------------------------------
## load("data/seu_obj_raw.RData")
## seu <- seu_obj


## -----------------------------------------------------------------------------
download.file("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz","10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")
untar("10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")
filt_mtx <- Seurat::Read10X("filtered_feature_bc_matrix/")


## ----eval=T, echo=F-----------------------------------------------------------
unlink("10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz", recursive=TRUE)
unlink("filtered_feature_bc_matrix/", recursive=TRUE)


## ----eval=T-------------------------------------------------------------------

filt_mtx_drop <- filt_mtx[rownames(filt_mtx) %in% names(amb),]
seu <- CreateSeuratObject(filt_mtx_drop)


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(filt_mtx_drop)


## ----eval=T-------------------------------------------------------------------
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
sce <- as.SingleCellExperiment(seu,assay = "RNA")


## -----------------------------------------------------------------------------
plotUMAP(sce, colour_by="seurat_clusters")


## -----------------------------------------------------------------------------
seu_list <- list()
seu_list[["woCorr"]] <- seu


## ----sec3_dropProc_ambProc_DropUtils,include=TRUE,eval=T----------------------
out <- removeAmbience(counts(sce), ambient=amb, groups=sce$seurat_clusters) 
rownames(out) <- rownames(sce)
colnames(out) <- colnames(sce)



## -----------------------------------------------------------------------------
seu_list[["withCorr"]] <- CreateSeuratObject(out)
seu_list[["withCorr"]] <- data_proc(seu_list[["withCorr"]])
seu_list[["withCorr"]] <- ScaleData(seu_list[["withCorr"]])
seu_list[["withCorr"]] <- quick_clust(seu_list[["withCorr"]])


## ----sec3_dropProc_testMarker_DropletUtil,include=TRUE,eval=T-----------------

mark_gene <- c("CCR7","CD8A","MS4A1","CD14","FCGR3A","FCER1A","GNLY","NKG7","PPBP")
mark_gene



## -----------------------------------------------------------------------------
DimPlot(seu_list$woCorr,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
DimPlot(seu_list$withCorr,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()



## -----------------------------------------------------------------------------
VlnPlot(seu_list$woCorr,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## -----------------------------------------------------------------------------
VlnPlot(seu_list$withCorr,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu_list, out, sce)
gc()


## ----sec3_dropProc_SoupX_prep,include=TRUE,eval=T-----------------------------

clust <- setNames(seu$seurat_clusters, Cells(seu))
seu_filt <- seu


## -----------------------------------------------------------------------------
library(SoupX)
soupOBJ <- SoupChannel(raw_mtx, filt_mtx)
soupOBJ <- setClusters(soupOBJ,clust)


## ----sec3_dropProc_SoupX_autoCor,include=TRUE,eval=F--------------------------
## soupOBJ <- autoEstCont(soupOBJ)
## 
## autoCor_mtx <- adjustCounts(soupOBJ)
## 


## ----eval=F, echo=F-----------------------------------------------------------
## save(soupOBJ, file="../data/soup.RData")
## save(autoCor_mtx, file="../data/auto_soup.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/soup.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/auto_soup.RData")


## -----------------------------------------------------------------------------
seu <- CreateSeuratObject(autoCor_mtx)
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
seu_autoCorr <- seu


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(seu)
gc()


## ----sec3_dropProc_SoupX_estCor,include=TRUE,eval=T---------------------------

mark_list <- list("CD4 T-cell"=c("IL7R","CCR7","S100A4"),"CD8 T-cell"=c("CD8A"),"B-Cell"=c("MS4A1"),
                  "Monocyte"=c("CD14","FCGR3A"),"DC"=c("FCER1A"),"NK"=c("NKG7","GNLY"),"Platelet"=c("PPBP"))

use_toEst <- estimateNonExpressingCells(soupOBJ, nonExpressedGeneList = mark_list)
soupOBJ <- calculateContaminationFraction(soupOBJ, mark_list, useToEst = use_toEst)
rho <- unique(soupOBJ$metaData$rho)
rho


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(use_toEst)
gc()


## ----eval=F-------------------------------------------------------------------
## soupOBJ <- setContaminationFraction(soupOBJ,rho,forceAccept=TRUE)
## estCor_mtx <- adjustCounts(soupOBJ)
## 


## ----eval=F, echo=F-----------------------------------------------------------
## save(estCor_mtx, file="../data/estCor_mtx.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/estCor_mtx.RData")


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(soupOBJ, rho)
gc()


## -----------------------------------------------------------------------------
seu <- CreateSeuratObject(estCor_mtx)
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
seu_estCorr <- seu


## ----sec3_dorpProc_SoupX_estCor_eval,include=TRUE,eval=T----------------------
mark_gene <- c("CCR7","CD8A","MS4A1","CD14","FCGR3A","FCER1A","GNLY","NKG7","PPBP")
mark_gene


## -----------------------------------------------------------------------------
DimPlot(seu_filt,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
VlnPlot(seu_filt, features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## -----------------------------------------------------------------------------

DimPlot(seu_autoCorr,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
VlnPlot(seu_autoCorr,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## -----------------------------------------------------------------------------
DimPlot(seu_estCorr,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
VlnPlot(seu_estCorr,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------
rm(seu_estCorr, seu_autoCorr)
gc()


## input_h5=the_raw_matrix_in_h5_format_from_cellranger #essential
## output_h5=assign_the_h5_file_path_for_the_cellbender_corrected_matrix # essential
## expect_cell=expected_cell_number_can_be_find_in_cellranger_Web_Summary # essential
## droplet_num=the_total_number_of_droplets_assigned_while_sequencing # default 25,000
## fpr=threshols_of_FALSE_POSITIVE_RATE # default 0.01
## epochs=number_of_epochs_to_train # default 150
## num_train=number_of_times_to_attempt_to_train_the_model # default 1. would speed up while setting greater
## #
## cellbender remove-background --input $input_h5 --output $output_h5 --expected-cells $expect_cell --total-droplets-included $droplet_num --fpr $fpr --epochs $epochs --num-training-tries $num_train --cuda False

## ----sec3_dropProc_cbFilt_loadData,include=TRUE,eval=T------------------------
cbFilt_mtx <- Read10X_h5("data/cbFilt_PBMCv3_20230324_filtered.h5")
dim(cbFilt_mtx)

dim(filt_mtx)


## ----sec3_dropProc_cb_Filt_dataPRoc,include=TRUE,eval=T-----------------------
message("processing matrix from CellRanger")
seu <- CreateSeuratObject(filt_mtx)
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
seu_filt <- seu

message("processing matrix from CellBender")
seu <- CreateSeuratObject(cbFilt_mtx)
seu <- data_proc(seu)
seu <- ScaleData(seu)
seu <- quick_clust(seu)
seu_cbFilt <- seu


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu)
gc()
  


## ----sec3_dorpProc_cbFilt_eval,include=TRUE,eval=T----------------------------
mark_gene <- c("CD3E","CCR7","CD8A","MS4A1","CD14","FCGR3A","FCER1A","GNLY","PPBP")
mark_gene


## -----------------------------------------------------------------------------
DimPlot(seu_filt,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
DimPlot(seu_cbFilt,group.by = "seurat_clusters",pt.size = 0.1,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------
VlnPlot(seu_filt,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)



## -----------------------------------------------------------------------------
VlnPlot(seu_cbFilt,features = mark_gene,group.by = "seurat_clusters",pt.size = 0)


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(bcrank, uniq, e.out, amb, filt_mtx, seu, sce, seu_list, out, soupOBJ, clust, autoCor_mtx, seu_autoCorr, use_toEst, rho, estCor_mtx, seu_estCorr, seu_cbFilt, seu_filt)
gc()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Trajectory and Psuedotime analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Trajectory and Psuedotime analysis

---
"    
  )
  
}



## ----sec3_trajAna_prep,include=TRUE,eval=T------------------------------------
sce <- readRDS("data/sceOBJ_Nestorowa_HSC.rds")
sce


## ----sec3_trajAna_buildTRAJ,include=TRUE,eval=T-------------------------------
library(slingshot)
sce.sling <- slingshot(sce,
                       cluster=sce$label,
                       reducedDim='PCA',
                       approx_points=100,
                       omega=TRUE)


## -----------------------------------------------------------------------------
colData(sce.sling)[1,1:5]


## ----sec3_trajAna_extTRAJ,include=TRUE,eval=T---------------------------------
embedded <- embedCurves(sce.sling, "UMAP")

embedded@metadata$lineages


## -----------------------------------------------------------------------------
pseudo.paths <- slingPseudotime(sce.sling)
head(pseudo.paths,2)



## -----------------------------------------------------------------------------
avg_pseudoTime <- rowMeans(pseudo.paths, na.rm=TRUE)
colData(sce.sling)$avg_pseudoTime <- avg_pseudoTime


## -----------------------------------------------------------------------------

plotUMAP(sce.sling,colour_by="avg_pseudoTime")


## ----sec3_trajAna_plotPR,include=TRUE,eval=T----------------------------------

embedded_curve <- slingCurves(embedded)

curve <- lapply(embedded_curve,function(x){
  dat <- data.frame(x$s[x$ord,]) # UMAP axis and the order of cells
  return(dat)})
names(curve) <- c("curve1","curve2","curve3")

head(curve$curve1)


## -----------------------------------------------------------------------------
plotUMAP(sce.sling,colour_by="avg_pseudoTime")+
  geom_path(data=curve$curve1,aes(x=UMAP1,y=UMAP2), color="blue")+
  geom_path(data=curve$curve2,aes(x=UMAP1,y=UMAP2),  color="black")+
  geom_path(data=curve$curve3,aes(x=UMAP1,y=UMAP2),  color="red")


## -----------------------------------------------------------------------------
embedded@metadata$lineages


## ----sec3_trajAna_refineTRAJ,include=TRUE,eval=T------------------------------
sce2 <- sce[,colData(sce)$label != 7]
sce.sling2 <- slingshot(sce2,cluster=sce2$label,reducedDim='PCA',approx_points=100,omega=TRUE)
pseudo.paths <- slingPseudotime(sce.sling2)
avg_pseudoTime <- rowMeans(pseudo.paths, na.rm=TRUE)
colData(sce.sling2)$avg_pseudoTime <- avg_pseudoTime


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(sce, sce2)
gc()


## -----------------------------------------------------------------------------

embedded <- embedCurves(sce.sling2, "UMAP")
embedded_curve <- slingCurves(embedded)
curve <- lapply(embedded_curve,function(x){
  dat <- data.frame(x$s[x$ord,]) # UMAP axis and the order of cells
  return(dat)})
names(curve) <- c("curve1","curve2")


## -----------------------------------------------------------------------------
plotUMAP(sce.sling2, colour_by="avg_pseudoTime") +
  geom_path(data=curve$curve1,aes(x=UMAP1,y=UMAP2), color="blue") +
  geom_path(data=curve$curve2,aes(x=UMAP1,y=UMAP2),  color="black")


## ----sec3_trajANA_trajSep,include=TRUE,eval=T---------------------------------
plotUMAP(sce.sling2, colour_by="label")


## -----------------------------------------------------------------------------
embedded@metadata$lineages$Lineage1
plotUMAP(sce.sling2,colour_by="slingPseudotime_1")+
  geom_path(data=curve$curve1,aes(x=UMAP1,y=UMAP2))


## -----------------------------------------------------------------------------
embedded@metadata$lineages$Lineage2
plotUMAP(sce.sling2, colour_by="slingPseudotime_2")+
  geom_path(data=curve$curve2,aes(x=UMAP1,y=UMAP2))


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(embedded, embedded_curve, curve)
gc()


## ----sec3_trajANA_driveL1, include=TRUE,eval=T--------------------------------
library(TSCAN)
res <- testPseudotime(sce.sling2, sce.sling2$slingPseudotime_1)


## -----------------------------------------------------------------------------
res$SYMBOL <- rownames(res)
res <- res[order(-res$logFC),]
head(res)



## -----------------------------------------------------------------------------
res$stat <- NA
res$stat[res$logFC > 0 & res$FDR < 0.05] <- "late"
res$stat[res$logFC < 0 & res$FDR < 0.05] <- "early"
res$stat[is.na(res$stat)] <- "nc"
res[1:2,]


## -----------------------------------------------------------------------------
table(res$stat)



## ----sec3_trajANA_driveL1_topEarly,eval=T-------------------------------------

sub <- res[res$stat=="early",]
sub <- sub[order(sub$logFC),]
top <- head(sub$SYMBOL,5)
top


## -----------------------------------------------------------------------------
plotExpression(sce.sling2,
               features = top, swap_rownames = "GENEID",
               x="slingPseudotime_1", colour_by = "label")


## ----sec3_trajANA_driveL1_topLate2,include=TRUE,eval=T------------------------

sub <- res[res$stat=="late",]
sub <- sub[order(-sub$logFC),]
top <- head(sub$SYMBOL,5)
plotExpression(sce.sling2,
               features = top,swap_rownames = "GENEID",
               x="slingPseudotime_1",colour_by = "label")


## ----sec3_trajANA_driveL1_topLate,include=TRUE,eval=T, echo=F-----------------

sub <- res[res$stat=="late",]
sub <- sub[order(-sub$logFC),]
top <- head(sub$SYMBOL,5)
plotExpression(sce.sling2,
               features = top,swap_rownames = "GENEID",
               x="slingPseudotime_1",colour_by = "label")


## ----sec3_trajANA_driveL2,include=TRUE,eval=T---------------------------------

res <- testPseudotime(sce.sling2,sce.sling2$slingPseudotime_2)
res$SYMBOL <- rownames(res)
res <- res[order(-res$logFC),]
res$stat <- NA
res$stat[res$logFC > 0 & res$FDR < 0.05] <- "late"
res$stat[res$logFC < 0 & res$FDR < 0.05] <- "early"
res$stat[is.na(res$stat)] <- "nc"

table(res$stat)


## ----sec3_trajANA_driveL2_topEarly,include=TRUE,eval=T------------------------
sub <- res[res$stat=="early",]
sub <- sub[order(sub$logFC),]
top <- head(sub$SYMBOL,5)
plotExpression(sce.sling2,
               features = top,swap_rownames = "GENEID",
               x="slingPseudotime_1",colour_by = "label")


## ----sec3_trajANA_driveL2_topLate,include=TRUE,eval=T-------------------------
sub <- res[res$stat=="late",]
sub <- sub[order(-sub$logFC),]
top <- head(sub$SYMBOL,5)
plotExpression(sce.sling2,
               features = top,swap_rownames = "GENEID",
               x="slingPseudotime_1",colour_by = "label")


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(sce.sling2)
rm(res, embedded, embedded_curve,sce2, pseudo.paths, avg_pseudoTime, curve, sce.sling, sce, tbl, panc_list, cell_type, mat, pred_res,sce_list,pred_cellType,anchors,ref_panc,feats,ref_list,seu_list,panc8, seu_obj,mnn_res,cell_type,snn.gr,cluster_mnn,universe,sce_list,sce,dec_list,dec,combined_dec,chosen_hvgs,hvgc_list,anchors,ifnb_merge,ifnb_list,feats,ifnb_list_rpca ,ifnb)
gc()




## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# CITE-Seq 

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# CITE-Seq 

---
"    
  )
  
}



## ----sec3_CITE_prep,include=TRUE,eval=T---------------------------------------

rna.mat <- readRDS("data/pbmc_umi_mtx.rds")
dim(rna.mat)

hto.mat <- readRDS("data/pbmc_hto_mtx.rds")
dim(hto.mat)
rownames(hto.mat)


## -----------------------------------------------------------------------------
seu_obj <- CreateSeuratObject(counts=rna.mat,project="citeSeq_demo")
seu_obj[["HTO"]] <- CreateAssayObject(counts=hto.mat)
seu_obj


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(hto.mat, rna.mat )
gc()


## ----sec3_CITE_clust,include=TRUE,eval=T--------------------------------------

DefaultAssay(seu_obj) <- "RNA"
seu_obj <- data_proc(seu_obj)
seu_obj <- ScaleData(seu_obj)


## -----------------------------------------------------------------------------
seu_obj <- quick_clust(seu_obj)
DimPlot(seu_obj,group.by = "seurat_clusters",pt.size = 0.2,label = TRUE)+NoLegend()


## ----sec3_CITE_hto,include=TRUE,eval=T----------------------------------------

seu_obj <- NormalizeData(seu_obj, assay="HTO", normalization.method="CLR")
seu_obj <- HTODemux(seu_obj, assay = "HTO", positive.quantile = 0.99)

head(seu_obj,2)


## -----------------------------------------------------------------------------
table(seu_obj$HTO_classification.global)
table(seu_obj$hash.ID)


## ----sec3_CITE_ridget,include=TRUE,eval=T-------------------------------------
RidgePlot(seu_obj,group.by = "hash.ID", assay = "HTO",
          features = rownames(seu_obj[["HTO"]])[1:2],ncol = 2)


## ----include=TRUE,eval=T------------------------------------------------------
RidgePlot(seu_obj,group.by = "hash.ID", assay = "HTO",
          features = rownames(seu_obj[["HTO"]])[3:4],ncol = 2)


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu_obj )
gc()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# SMART-Seq 

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# SMART-Seq 

---
"    
  )
  
}



## ----SMARTSeq_loadMTX,include=TRUE,eval=TRUE----------------------------------
mtx <- read.delim("data/GSE151334_counts.mouse.tsv.gz")
mtx[1:2,1:5]


## ----SMARTSeq_loadSeu,include=TRUE,eval=TRUE----------------------------------
library(Seurat)
seu <- CreateSeuratObject(mtx,project="GSE151334",min.cells=3,min.features=300)
seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern = "^mt-")
seu <- SCTransform(seu,variable.features.n = 2000)



## ----echo=F, eval=F-----------------------------------------------------------
## 
## # seu[["cellType"]] <- gsub("(.*)_(.*)_(.*)_(.*)_(.*)","\\1_\\2",Cells(seu))
## # seu[["date"]] <- gsub("(.*)_(.*)_(.*)_(.*)_(.*)","\\3",Cells(seu))


## ----echo=F, eval=F-----------------------------------------------------------
## library(reticulate)
## # py_path <- "//Users/mattpaul/myMiniconda/bin/python3"
## # Sys.setenv(RETICULATE_PYTHON=py_path)


## ----SMARTSeq_doublet,include=TRUE,eval=F-------------------------------------
## library(reticulate)
## reticulate::py_install("scrublet")
## scr <- reticulate::import("scrublet")
## 
## mat <- GetAssayData(seu,assay = "RNA",slot = "counts")
## mat <- as.matrix(mat)
## scr <- reticulate::import("scrublet")
## scrub <- scr$Scrublet(t(mat))
## doublet <- scrub$scrub_doublets()
## names(doublet) <- c("doublet_score","doublet")
## seu[["doublet_score"]] <- doublet$doublet_score
## seu[["doublet"]] <- doublet$doublet


## ----SMARTSeq_estCC,include=TRUE,eval=F,echo=FALSE----------------------------
## 
## cc_ge <- read.delim("data/mouse_regev_lab_cell_cycle_genes.txt",header = FALSE,stringsAsFactors = FALSE)
## ccGene_list <- list("S"=cc_ge$V1[1:43],"G2M"=cc_ge$V1[44:length(cc_ge)])
## seu <- CellCycleScoring(seu,ccGene_list$S,ccGene_list$G2M)


## ----SMARTSeq_norm,include=TRUE,eval=F----------------------------------------
## set.seed(42)
## seu <- RunPCA(seu,npcs = 50,verbose = FALSE)
## seu <- FindNeighbors(seu,dims = 1:15,reduction = "pca")
## seu <- RunUMAP(seu,dims = 1:15,reduction = "pca")
## seu <- FindClusters(seu,resolution = 0.5)
## 


## ----echo=F, eval=F-----------------------------------------------------------
## saveRDS(seu,"data/seuOBJ_mouseEB_ori.rds"")
## 


## ----echo=F, eval=T-----------------------------------------------------------
seu <- readRDS("data/seuOBJ_mouseEB_ori.rds")



## ----SMARTSeq_quickClust,include=TRUE,eval=TRUE-------------------------------
library(ggplot2)
VlnPlot(seu,
        features = c("nCount_RNA","nFeature_RNA","percent.mt","doublet_score","S.Score","G2M.Score"),
        group.by = "seurat_clusters")


## -----------------------------------------------------------------------------
FeatureScatter(seu,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "seurat_clusters")


## -----------------------------------------------------------------------------
FeatureScatter(seu,feature1 = "nCount_RNA",feature2 = "percent.mt",group.by = "seurat_clusters") + geom_hline(yintercept = 10,color="black",linetype="dashed")


## -----------------------------------------------------------------------------
table(seu$doublet=="TRUE" | seu$percent.mt > 10)


## ----SMARTSeq_filt_quickClust,include=TRUE,eval=TRUE--------------------------
set.seed(42)
seu <- subset(seu,subset=percent.mt <= 10 & doublet == "FALSE")
seu <- SCTransform(seu,variable.features.n = 2000,vars.to.regress = c("percent.mt"))
seu <- RunPCA(seu,npcs = 50,verbose = FALSE)
seu <- RunUMAP(seu,dims = 1:15,reduction = "pca")
seu <- FindNeighbors(seu,dims = 1:15,reduction = "pca")
seu <- FindClusters(seu,resolution = 0.5)



## -----------------------------------------------------------------------------
DimPlot(seu,group.by = "seurat_clusters",pt.size = 0.2,label = TRUE)+NoLegend()


## ----SMARTseq_findMark,include=TRUE,eval=F------------------------------------
## markers <- FindAllMarkers(seu,only.pos = TRUE,min.pct = 0.25,
##                           logfc.threshold = 0.25)
## library(dplyr)
## top_genes <- markers %>% group_by(cluster) %>%
##   slice_max(n=5,order_by=avg_log2FC)


## ----eval=F, echo=F-----------------------------------------------------------
## 
## write.csv(top_genes,"data/top_genes.csv")


## ----echo=F-------------------------------------------------------------------
top_genes <- read.csv("data/top_genes.csv")


## -----------------------------------------------------------------------------

DoHeatmap(seu,features=c(top_genes$gene))


## ----SMARTSeq_markExp,include=TRUE,eval=TRUE----------------------------------
VlnPlot(seu,features = c("Nanog","Pou5f1"),
        ncol = 2,
        pt.size = 0)


## ----SMARTSeq_markExp2,include=TRUE,eval=TRUE---------------------------------
VlnPlot(seu,features = c("Pax6","Stau2"),
        ncol = 2,
        pt.size = 0)


## ----SMARTSeq_markExp3,include=TRUE,eval=TRUE---------------------------------
VlnPlot(seu,features = c("Epcam","Gata6"),
         ncol = 2,
         pt.size = 0)


## ----SMARTSeq_markExp4,include=TRUE,eval=TRUE---------------------------------
VlnPlot(seu,features = c("Dlk1","Postn"),
         ncol = 2,
         pt.size = 0)


## -----------------------------------------------------------------------------
seu[["layers"]] <- NA
seu$layers[seu$seurat_clusters %in% c(0,4,5)] <- "primary_ES"
seu$layers[seu$seurat_clusters %in% c(2,3,10,11,12)] <- "ectoderm"
seu$layers[seu$seurat_clusters %in% c(1,9)] <- "mesoderm"
seu$layers[seu$seurat_clusters %in% c(6,7,8)] <- "endoderm"


## -----------------------------------------------------------------------------
DimPlot(seu,group.by = "layers",pt.size = 0.2)

