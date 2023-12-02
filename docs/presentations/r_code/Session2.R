params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(bioMart))
suppressPackageStartupMessages(require(SeuratWrappers))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----warning=F, message=F, echo=F---------------------------------------------
library(RCurl)
if(!url.exists("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")){stop("Download path broken")}



## ----eval=F-------------------------------------------------------------------
## download.file("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz","10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")


## ----eval=F-------------------------------------------------------------------
## untar("10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")
## 


## ----eval=F-------------------------------------------------------------------
## dir()
## 


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



## ----include=F,echo=T,eval=TRUE-----------------------------------------------
mtx_dir <- "filtered_feature_bc_matrix"



## ----load_data,include=TRUE,eval=F--------------------------------------------
## library(Seurat)
## 
## mtx <- Seurat::Read10X(mtx_dir)
## 


## ----eval=F, echo=F-----------------------------------------------------------
## a <- is(mtx)
## b <- head(mtx)
## save(a,b, file="data/mtx_res.RData")


## ----eval=F, echo=F-----------------------------------------------------------
## library(Seurat)
## load("data/seurat_read.RData")


## ----eval=T, echo=F-----------------------------------------------------------
library(Seurat)
load("data/mtx_res.RData")


## ----eval=F, echo=T-----------------------------------------------------------
## is(mtx)


## ----eval=T, echo=F-----------------------------------------------------------
a


## ----eval=F, echo=T-----------------------------------------------------------
## head(mtx)


## ----eval=T, echo=F-----------------------------------------------------------
b


## ----load_h5,include=TRUE,eval=FALSE------------------------------------------
## h5_file <- "path to matrix h5 file"
## h5_file <- "~/Downloads/10k_PBMC_3p_nextgem_Chromium_Controller_molecule_info.h5"
## mtx <- Seurat::Read10X_h5(h5_file)


## ----load_CreateOBJ,include=TRUE,eval=TRUE------------------------------------
sample_id <- "PBMC_10k" # sample name
min_gene <- 200
min_cell <- 10 


## ----eval=F-------------------------------------------------------------------
## seu_obj <- Seurat::CreateSeuratObject(mtx, project=sample_id, min.cells=min_cell, min.features=min_gene)


## ----eval=F, echo=F-----------------------------------------------------------
## save(seu_obj, file="data/seu_obj_raw.RData")


## ----eval=T, echo=F-----------------------------------------------------------
load("data/seu_obj_raw.RData")


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
mt_gene_det <- mt_gene[mt_gene %in% rownames(seu_obj)]
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, features = mt_gene_det)
summary(seu_obj$percent.mt)


## ----load_estMT,include=TRUE,eval=TRUE----------------------------------------
seu_obj[["percent.mt2"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
summary(seu_obj$percent.mt2)



## ----eval =F------------------------------------------------------------------
## library(ggplot2)
## dat <- data.frame(byPattern=seu_obj$percent.mt, byGene=seu_obj$percent.mt2,stringsAsFactors = FALSE)
## cor_val <- cor.test(dat$byPattern,dat$byGene,method = "spearman")
## ggplot(dat,aes(x=byPattern,y=byGene))+geom_point()+geom_smooth()+
##   labs(x="% of MT, est by genes",y="% of MT,est by pattern",
##        subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
##   theme_classic()


## ----echo =F, fig.height=4,fig.width=7----------------------------------------
library(ggplot2)
dat <- data.frame(byPattern=seu_obj$percent.mt, byGene=seu_obj$percent.mt2,stringsAsFactors = FALSE)
cor_val <- cor.test(dat$byPattern,dat$byGene,method = "spearman")
ggplot(dat,aes(x=byPattern,y=byGene))+geom_point()+geom_smooth()+
  labs(x="% of MT, est by genes",y="% of MT,est by pattern",
       subtitle = paste0("rho=",round(cor_val$estimate,3),"; p-value=",cor_val$p.value[1]))+
  theme_classic()


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
## top10 <- head(VariableFeatures(seu_obj), 10)
## 
## plot1 <- VariableFeaturePlot(seu_obj)
## plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## plot2


## ----echo=F,fig.height=4,fig.width=7------------------------------------------
top10 <- head(VariableFeatures(seu_obj), 10)

plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


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
DefaultAssay(seu_obj) <- "RNA"
seu_obj <- CellCycleScoring(seu_obj, s.features = feat_s, g2m.features = feat_g2m)



## -----------------------------------------------------------------------------

dat_s <- data.frame(cell_id=Cells(seu_obj), cat="S_Score", Phase=seu_obj$Phase, score=seu_obj$S.Score)

dat_g2m <- data.frame(cell_id=Cells(seu_obj), cat="G2M_Score", Phase=seu_obj$Phase, score=seu_obj$G2M.Score)

dat <- rbind(dat_s, dat_g2m)
dat$Phase <- factor(dat$Phase, levels = c("G1","S","G2M"))


## ----fig.height=4,fig.width=7-------------------------------------------------

ggplot(dat,aes(x=Phase, y=score, fill=Phase))+geom_boxplot()+
  labs(x="",y="Score",fill="Phase")+
  facet_wrap(~cat)+theme_classic()


## ----ccPhase_cyclon_prep,include=TRUE,eval=TRUE-------------------------------
library(scran)

sce <- as.SingleCellExperiment(seu_obj,assay = "RNA")
rowData(sce)$SYMBOL <- rownames(sce)
sce


## -----------------------------------------------------------------------------

load("data/ccGene_mouse_human_geneSymbol_ensemblID_20220817.RData")
ccGene_hs <- ccGene_mm_hs$human_symbol
lapply(ccGene_hs, function(x){head(x,2)})


## ----ccPhase_cyclon_proc,include=TRUE,eval=FALSE------------------------------
## assignments <- cyclone(sce, ccGene_hs, gene.names=rowData(sce)$SYMBOL)
## 


## ----echo=F-------------------------------------------------------------------
#save(assignments, file="data/cyclone.RData")


## ----echo=T-------------------------------------------------------------------

load("data/cyclone.RData")


## -----------------------------------------------------------------------------
lapply(assignments, head)


## -----------------------------------------------------------------------------
seu_obj[["cyclon_Phase"]] <- assignments$phases
seu_obj[["cyclon_G1Score"]] <- assignments$scores$G1
seu_obj[["cyclon_SScore"]] <- assignments$scores$S
seu_obj[["cyclon_G2MScore"]] <- assignments$scores$G2M


## ----ccPhase_cyclon_boxPlot,include=TRUE,eval=T-------------------------------
dat_g1 <- data.frame(cell_id=Cells(seu_obj), cat="cyclon_G1Score", Phase=seu_obj$cyclon_Phase, score=seu_obj$cyclon_G1Score)

dat_s <- data.frame(cell_id=Cells(seu_obj), cat="cyclon_SScore", Phase=seu_obj$cyclon_Phase, score=seu_obj$cyclon_SScore)

dat_g2m <- data.frame(cell_id=Cells(seu_obj), cat="cyclon_G2MScore", Phase=seu_obj$cyclon_Phase, score=seu_obj$cyclon_G2MScore)

dat <- rbind(dat_g1,dat_s,dat_g2m)

dat$Phase <- factor(dat$Phase,levels = c("G1","S","G2M"))
dat$cat <- factor(dat$cat,levels = c("cyclon_G1Score","cyclon_SScore","cyclon_G2MScore"))


## ----fig.height=4,fig.width=7-------------------------------------------------
ggplot(dat,aes(x=Phase,y=score,fill=Phase))+geom_boxplot()+
  labs(x="",y="Score",fill="Phase")+
  facet_wrap(~cat)+theme_classic()


## -----------------------------------------------------------------------------
table(seu_obj$Phase)


## -----------------------------------------------------------------------------
table(seu_obj$cyclon_Phase)


## -----------------------------------------------------------------------------
table(seu_obj$Phase,seu_obj$cyclon_Phase)


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
## # library(Herper)
## #
## # conda_install  <- install_CondaTools("scrublet", "scRNA", pathToMiniConda = "../mini")
## #
## # Sys.setenv('RETICULATE_PYTHON'=file.path(conda_install$pathToEnvBin, "python"))


## ----eval=F-------------------------------------------------------------------
## library(reticulate)
## reticulate::py_install("scrublet")
## 
## scr <- reticulate::import("scrublet")
## 


## ----det_doublet_est,include=TRUE,eval=F--------------------------------------
## 
## mat <- GetAssayData(seu_obj, assay = "RNA", slot = "counts")
## mat <- as.matrix(mat)
## 


## ----eval=F-------------------------------------------------------------------
## # Loading the scrublet library
## scr <- import("scrublet")
## # Run scrublet
## scrub <- scr$Scrublet(t(mat))
## # Extract scrublet results
## doublet <- scrub$scrub_doublets()
## names(doublet) <- c("doublet_score","doublet")


## ----eval=F, echo=F-----------------------------------------------------------
## #save(doublet,file = "data/doublet.RData")
## 


## ----eval=T, echo=F-----------------------------------------------------------

load("data/doublet.RData")



## -----------------------------------------------------------------------------
summary(doublet$doublet_score)


## -----------------------------------------------------------------------------
table(doublet$doublet)


## ----det_doublet_pres,include=TRUE,eval=TRUE----------------------------------
seu_obj[["doublet_score"]] <- doublet$doublet_score
seu_obj[["doublet"]] <- doublet$doublet


## ----fig.height=4,fig.width=7-------------------------------------------------

VlnPlot(seu_obj, group.by = "doublet",
        features = c("nCount_RNA","nFeature_RNA","doublet_score"),
        pt.size = 0)


## ----fig.height=4,fig.width=7-------------------------------------------------
FeatureScatter(seu_obj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",pt.size = 0.1,group.by = "doublet")


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
## VlnPlot(seu_obj, group.by = "dset",
##         features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
##         pt.size = 0)


## ----include=TRUE,echo=F, fig.height=4,fig.width=7----------------------------
VlnPlot(seu_obj, group.by = "dset",
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        pt.size = 0)


## ----qcPlot_scatter_pres1,include=TRUE,eval=F---------------------------------
## 
## FeatureScatter(seu_obj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",
##                group.by = "doublet")


## ----include=TRUE,eval=T, echo=F, fig.height=4,fig.width=7--------------------

FeatureScatter(seu_obj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",
               group.by = "doublet")


## ----include=TRUE,eval=F------------------------------------------------------
## FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")


## ----qcPlot_scatter_pres2,include=TRUE,eval=TRUE, echo=F, fig.height=4,fig.width=7----
FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")


## ----include=TRUE,eval=F------------------------------------------------------
## RidgePlot(seu_obj,group.by = "doublet",features = c("doublet_score"))


## ----qcPlot_ridgePlot_pres,include=TRUE,eval=TRUE, echo=F, fig.height=4,fig.width=7----
RidgePlot(seu_obj,group.by = "doublet",features = c("doublet_score"))


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
table(seu_obj$doublet=="TRUE" | seu_obj$percent.mt >= 10)


## ----eval=F-------------------------------------------------------------------
## seu_filt <- subset(seu_obj, subset=doublet=="FALSE" &
##                      percent.mt < 10)
## 


## ----echo=F, eval=TRUE--------------------------------------------------------
#save(seu_filt,file="data/seu_filt.RData")
rm(doublet)
rm(seu_obj)
#load("data/seu_filt.RData")


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
## pot_conf <- c("percent.mt","doublet_score","cyclon_G1Score","cyclon_SScore","cyclon_G2MScore","cyclon_Phase")
## seu_filt <- ScaleData(seu_filt, vars.to.regress = pot_conf)
## 


## ----eval=F, echo=F-----------------------------------------------------------
## #save(seu_filt, sce, file="data/SCT2.RData")
## save(seu_filt, file="data/SCT2.RData")
## #saveRDS(seu_filt, file="data/SCT2.rds")
## #save(sce, file="data/sce.RData")
## 


## ----sceload, eval=T, echo=F--------------------------------------------------
#load("data/sce.RData")


## ----sctload, eval=T, echo=F--------------------------------------------------
load("data/SCT2.RData")
#seu_filt <- load("data/SCT2.rds")


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
DefaultAssay(seu_filt) <- "RNA"
seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)



## ----dimRed_elbow,include=TRUE,eval=TRUE, fig.height=4,fig.width=7------------
ElbowPlot(seu_filt, ndims = 50, reduction = "pca")


## ----dimRed_elbow2,include=TRUE,eval=T----------------------------------------
library(dplyr)
# Determine percent of variation associated with each PC
pct <- seu_filt[["pca"]]@stdev / sum(seu_filt[["pca"]]@stdev) * 100
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
## plot_df <- data.frame(pct = pct, cumu = cumu,rank = 1:length(pct))
## 
## ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
##   geom_vline(xintercept = 90, color = "grey") +
##   geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
##   theme_bw()


## ----eval=T, echo=F, fig.height=4,fig.width=7---------------------------------
plot_df <- data.frame(pct = pct, cumu = cumu,rank = 1:length(pct))

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


## ----dimRed_UMAP,include=TRUE,eval=T------------------------------------------
seu_filt <- FindNeighbors(seu_filt,dims = 1:pc,reduction = "pca")

seu_filt <- RunUMAP(seu_filt,dims = 1:pc,reduction = "pca")


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt)



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



## ----clust_default,include=TRUE,eval=T----------------------------------------
seu_filt <- FindClusters(seu_filt, resolution = 0.5)
seu_filt[["cluster_byDefault"]] <- seu_filt$seurat_clusters


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2)+NoLegend()




## ----clust_resoEval,include=TRUE,eval=T---------------------------------------
library(clustree)

reso <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
reso_res <- lapply(1:length(reso), function(x,seu_filt,reso){
  seu_filt <- FindClusters(seu_filt,resolution = reso[x])
  clust <- setNames(seu_filt$seurat_clusters,Cells(seu_filt))
  return(clust)}, seu_filt, reso)
names(reso_res) <- paste0("k",1:length(reso))


## ----eval=F-------------------------------------------------------------------
## remotes::install_github("thomasp85/tweenr")


## ----include=TRUE,eval=T------------------------------------------------------
k_tab <- do.call(cbind,reso_res)
k_dat <- as.data.frame(k_tab)

head(k_dat,2)


## ----clust_resoEval2,include=TRUE,eval=F--------------------------------------
## clustree(k_dat, prefix = "k", node_colour = "sc3_stability")


## ----include=TRUE,eval=T, echo=F, warning=FALSE, fig.height=4,fig.width=7-----
clustree(k_dat, prefix = "k", node_colour = "sc3_stability")


## ----clust_resoFix,include=TRUE,eval=T, fig.height=4,fig.width=7, warning=FALSE, message=FALSE----
seu_filt <- FindClusters(seu_filt, resolution = 0.4)


## ----fig.height=4,fig.width=7-------------------------------------------------
DimPlot(seu_filt, group.by = "seurat_clusters",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Optimized Clusters")


## ----include=TRUE,eval=T, fig.height=4,fig.width=7----------------------------
DimPlot(seu_filt, group.by = "cluster_byDefault",label = TRUE,pt.size = 0.2) + NoLegend() + ggtitle("Default Clusters")


## ----markGene_cal,include=TRUE,eval=FALSE-------------------------------------
## markers <- FindAllMarkers(seu_filt, only.pos = TRUE,
##                           min.pct = 0.25, logfc.threshold = 0.25)


## ----eval=T, echo=F-----------------------------------------------------------
#save(markers, file="data/markers.RData")
load("data/markers.RData")


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
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(5)] <- "B-cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(1,4,6)] <- "Naive CD4 T-cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(3)] <- "Memory CD4 T-cells"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(7)] <- "CD8 T-cells/NKs"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(8)] <- "NKs"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(2)] <- "FCGR3A+ Monocytes"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(0)] <- "CD14+ Monocytes"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(9)] <- "Platelets"
seu_filt$cellType_byClust[seu_filt$seurat_clusters %in% c(10)] <- "DCs"
seu_filt$cellType_byClust[is.na(seu_filt$cellType_byClust)] <- "Unspecified"


## -----------------------------------------------------------------------------

table(seu_filt$cellType_byClust)


## ----fig.height=4,fig.width=7-------------------------------------------------

DimPlot(seu_filt, group.by = c("seurat_clusters","cellType_byClust"),label = TRUE,pt.size = 0.2)+NoLegend()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential gene expression between cell types

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html>

---
"
  )
}else{
  cat("# Differential gene expression between cell types


---
"
  )

}



## ----sec2_degCal,include=TRUE,eval=T------------------------------------------
deg <- FindMarkers(seu_filt, group.by = "cellType_byClust",
                   ident.1 = "CD14+ Monocytes",ident.2 = "FCGR3A+ Monocytes",
                   logfc.threshold = 0.25,
                   test.use = "wilcox",
                   min.pct = 0.1)




## -----------------------------------------------------------------------------
deg$Symbol <- rownames(deg)
deg <- deg[order(-deg$avg_log2FC),]
head(deg)



## ----sec2_degCal_volcano,include=TRUE,eval=T----------------------------------
deg$sig <- NA
deg$sig[deg$avg_log2FC >= 0.585 & deg$p_val_adj < 0.05] <- "CD14+"
deg$sig[deg$avg_log2FC <=- 0.585 & deg$p_val_adj < 0.05] <- "FCGR3A+"
deg$sig[is.na(deg$sig)] <- "nc"
deg$sig <- factor(deg$sig, levels = c("CD14+","nc","FCGR3A+"))

table(deg$sig)


## ----fig.height=4,fig.width=7-------------------------------------------------
ggplot(deg, aes(x=avg_log2FC, y=-log10(p_val_adj), color=sig)) + geom_point() +
  scale_color_manual(values = c("red", "grey", "blue")) +
  labs(x="log2_Fold-Changes", y="-log10_adjPV", color="") +
  theme_classic()


## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(seu_filt)
rm(markers)
rm(sce)
rm(deg)
rm(cellType_dat)
rm(umap)
rm(umap_tab)
rm(meta)
rm(clust_dat)

unlink("umap_dat.csv")
unlink("cellType_dat.csv")

rm(clust)
rm(mat)
rm(k_dat)
rm(plot_df)
rm(dat)
rm(assignments)
rm(ccGene_hs)
rm(sct_mat)
rm(log_mat)
rm(log_avgExp )
rm(sct_avgExp)
gc()

unlink("clust_dat.csv")


