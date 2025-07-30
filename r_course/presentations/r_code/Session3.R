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


## ----sec3_loadPack,include=F,eval=T-------------------------------------------
library(Seurat)
library(scran)
library(scater)
library(ggplot2)
library(pheatmap)
library(TSCAN)
library(SoupX)
library(DropletUtils)
library(scuttle)

library(ggpubr)
library(DESeq2)
library(tibble)
library(dplyr)
library(MAST)
library(tidyr)


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
pheatmap(tbl, scale = "row", main="Harmony", treeheight_row = 0, , treeheight_col = 0)


## ----echo=F-------------------------------------------------------------------
  tbl <- table(my_seu_merge_rpca$sample_id, my_seu_merge_rpca$seurat_clusters)
pheatmap(tbl, scale = "row", main="rPCA", treeheight_row = 0, , treeheight_col = 0)



## ----echo=F-------------------------------------------------------------------
tbl <- table(seu_merge_harmony$sample_id, seu_merge_harmony$seurat_clusters)
pheatmap(tbl, scale = "row", main="Harmony", treeheight_row = 0, , treeheight_col = 0)


## ----eval=T, echo=F-----------------------------------------------------------
rm(seu_merge )
rm(seu_merge_harmony)



## ----echo=F, warning=FALSE, message=F-----------------------------------------
rm(seu_merge_harmony, my_seu_list_rpca,my_seu_list,seu_merge)
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



## ----echo=F-------------------------------------------------------------------

DimPlot(my_seu_merge_rpca, group.by = "paper_annot")


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


## ----echo=F, warning=FALSE, message=F-----------------------------------------
rm(hpcad)
gc()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential expression analysis

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Differential expression analysis

---
"    
  )
  
}



## ----echo=T, eval=F-----------------------------------------------------------
# # downloaded from dropbox
# seu_obj <- readRDS("~/Downloads/integrated.rds")
# seu_obj


## ----echo=F, eval=T-----------------------------------------------------------
seu_obj <- readRDS("./data/integrated.rds")
seu_obj


## ----results='hide',include=T,echo=T, eval=T, fig.height=4--------------------

p_all <- DimPlot(seu_obj,  reduction = "umap", group.by = "paper_annot", label = T) 
p_all


## ----seu_setup,include=T,echo=T, fig.height=3, fig.width=6--------------------
# metadata column for condition
seu_obj$condition <- ifelse(grepl("AD", seu_obj$sample_id), "AD", "CON")

# metadata column for condition + cell type
seu_obj$celltype_condition <- paste(seu_obj$paper_annot, seu_obj$condition, sep = "_")

table(seu_obj$celltype_condition)



## -----------------------------------------------------------------------------

plot_percents <- data.frame(table(seu_obj$celltype_condition)) %>%
  tidyr::separate(Var1, sep = "_", into = c("cell_type", "condition")) %>%
  group_by(cell_type) %>%
  mutate(percent = Freq/sum(Freq))

ggplot(plot_percents, aes(x = cell_type, y = percent, fill = condition)) +
  geom_bar(stat = "identity") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----results='hide',include=T,echo=T, eval=T, fig.width=10, fig.height=5------

DimPlot(subset(seu_obj, paper_annot == "Excitatory Neurons"),  reduction = "umap", group.by = "condition", label = F) + ggtitle("Excitatory Neurons only - AD vs Control")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Prepare object

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Prepare object

---
"    
  )
  
}



## ----results='asis',include=T,echo=T, eval=T----------------------------------
DefaultAssay(seu_obj)
DefaultAssay(seu_obj) <- "RNA"
DefaultAssay(seu_obj)


## ----results='asis',include=T,echo=T, eval=T----------------------------------
Layers(seu_obj)


## ----results='asis',include=T,echo=T, eval=T----------------------------------
seu_obj <- JoinLayers(seu_obj)
Layers(seu_obj)


## ----setup3,include=T,echo=T, fig.height=3------------------------------------
counts <- GetAssayData(seu_obj, assay = "RNA", layer = "counts")
counts_mat <- as.matrix(counts)
actin_counts <- counts_mat[rownames(counts_mat) == "ACTB"]
hist(actin_counts, breaks = 50)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Wilcoxon Rank Sum Test

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Wilcoxon Rank Sum Test

---
"    
  )
  
}



## ----wilcox,include=T,echo=T--------------------------------------------------
library(tibble)
library(dplyr)
Idents(seu_obj) <- "celltype_condition"
de_wilcox <- FindMarkers(seu_obj, ident.1 = "Excitatory Neurons_AD", ident.2 = "Excitatory Neurons_CON", logfc.threshold = 0)  %>% 
  tibble::rownames_to_column("geneID")

head(de_wilcox, 3)


## ----wilcox_volcano,include=T,echo=T------------------------------------------
de_wilcox$sig <- de_wilcox$p_val_adj < 0.05

ggplot(de_wilcox, aes(x = avg_log2FC, y = -log10(p_val), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# DE using parametric methods

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# DE using parametric methods

---
"    
  )
  
}



## ----include=T,echo=T, fig.height=3-------------------------------------------

counts_seu <- GetAssayData(seu_obj, assay = "RNA", layer = "data")
# SNX31 was identified as differentially expressed in excitatory neurons
seu_obj$SNX31 <- counts_seu[rownames(counts_seu) == "SNX31", ]
exn_toPlot <- seu_obj@meta.data[seu_obj@meta.data$paper_annot == "Excitatory Neurons", ]
ggplot(exn_toPlot, aes(x = SNX31, fill = condition)) + geom_density() + ggtitle("SNX31 expression in AD vs CON in Excitatory Neruons") + theme_classic() + xlim(0,2)


## ----mast_prep, include=T,echo=T, eval=T--------------------------------------
library(MAST)
seu_exn <- subset(seu_obj, paper_annot == "Excitatory Neurons")
seu_exn_data <- GetAssayData(seu_exn, assay = "RNA", layer = "data")
# create SingleCellAssay for MAST
sca_exn <- MAST::FromMatrix(exprsArray = as.matrix(seu_exn_data),
                           cData = seu_exn@meta.data,
                           fData = data.frame(gene_id = rownames(seu_exn))
)
sca_exn


## ----mast_prep2, include=T,echo=T, eval=T-------------------------------------
cdr <- colSums(assay(sca_exn)>0)
cdr_scale <- scale(cdr)
head(cdr_scale)


## ----mast_prep3, include=T,echo=T, eval=T-------------------------------------

colData(sca_exn)$ngeneson <- cdr_scale
cond <- factor(colData(sca_exn)$condition,
               levels = c("AD","CON"))
cond <- relevel(cond,"CON")
colData(sca_exn)$Group <- cond

colData(sca_exn)[1:3, c("paper_annot", "ngeneson", "Group")]



## ----mast_cdr_pca-------------------------------------------------------------
pca_embed <- data.frame(Embeddings(seu_exn, reduction = "pca")) # get PCA embeddings

# divide variance from each PC by total variance to get amount of variance explained by that PC
pca_vars <-  Stdev(seu_exn, reduction = "pca")^2
total_var <- sum(pca_vars) 
pct_var_explained <- pca_vars/total_var * 100 
pct_var_explained <- round(pct_var_explained, 1)

for_plot <- data.frame(colData(sca_exn))
for_plot <- merge(for_plot, pca_embed, by = 0)


## ----mast_cdr_pca_plot,fig.width=7--------------------------------------------
library(ggpubr)
p1 <- ggplot(for_plot, aes(x = ngeneson, y = PC_1, color = condition)) + geom_point() + ylab(paste0("PC_1 (", pct_var_explained[1],"%)")) 
p2 <- ggplot(for_plot, aes(x = ngeneson, y = PC_2, color = condition)) + geom_point() + ylab(paste0("PC_2 (", pct_var_explained[2],"%)")) 
ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend="bottom")


## ----mast_noRE,include=T,echo=T, eval=F---------------------------------------
# 
# zlmCond_exn <- zlm(~ngeneson + Group, sca_exn)
# sumZLM_exn <- summary(zlmCond_exn,
#                   doLRT=paste0("Group","AD"))
# sumDT_exn <- sumZLM_exn$datatable


## ----mast_noRE_saveDT,include=T,echo=F, eval=F--------------------------------
# saveRDS(sumDT_exn, "sumDT_exn.rds")


## ----include=T,echo=T, eval=F-------------------------------------------------
# sumDT_exn <- readRDS("~/Downloads/sumDT_exn.rds")


## ----mast_noRE_loadDT,include=T,echo=F, eval=T--------------------------------
sumDT_exn <- readRDS("./data/sumDT_exn.rds")


## -----------------------------------------------------------------------------
sumDT_exn[1:12, c("primerid", "component", "contrast", "Pr(>Chisq)", "coef")]


## ----mast_noRE_getRes,include=T,echo=T, eval=T--------------------------------
de_mast_exn <- merge(sumDT_exn[sumDT_exn$component=="H"&sumDT_exn$contrast==paste0("Group","AD"),],
               sumDT_exn[sumDT_exn$component=="logFC"&sumDT_exn$contrast==paste0("Group","AD"),],by='primerid')
de_mast_exn <- dplyr::select(de_mast_exn,primerid,coef.y,z.y,`Pr(>Chisq).x`)
names(de_mast_exn) <- c("geneID","log2Fc","z","pvalue")
de_mast_exn$FDR <- p.adjust(de_mast_exn$pvalue,'fdr')
de_mast_exn <- de_mast_exn[order(de_mast_exn$FDR),]

head(de_mast_exn, 3)


## ----mast_noRE_plot,include=T,echo=T, eval=T----------------------------------
de_mast_exn$sig <- de_mast_exn$FDR < 0.05

ggplot(de_mast_exn, aes(x = log2Fc, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw()


## ----mast_topSig, echo=T, eval=T, fig.width=8---------------------------------
de_mast_exn_df <- data.frame(de_mast_exn)
de_mast_exn_df <- na.omit(de_mast_exn_df)
top_down_AD <- head(de_mast_exn_df[de_mast_exn_df$log2Fc < 0, ], 5)

VlnPlot(seu_exn, features = top_down_AD$geneID, stack = T, flip = T, pt.size = 1) + scale_x_discrete(labels=c('CON', 'AD')) + NoLegend()



## ----mast_sorta_down, echo=T, eval=T, fig.width=10----------------------------
de_mast_exn_sig <- de_mast_exn_df[de_mast_exn_df$FDR < 0.05, ]
sorta_down_AD <- tail(de_mast_exn_sig[de_mast_exn_sig$log2Fc < 0, ], 5)
sorta_down_AD



## ----mast_violin2, echo=T, eval=T, fig.width=7--------------------------------

VlnPlot(seu_exn, features = sorta_down_AD$geneID, stack = T, flip = T, pt.size = 1) + scale_x_discrete(labels=c('CON', 'AD')) + NoLegend()




## ----mast_xist_feature, echo=T, eval=T, fig.width=10, fig.height=3------------
FeaturePlot(seu_exn, features = "XIST", split.by = "sample_id")


## ----mast_sex_cov1, eval=T, echo=T--------------------------------------------

colData(sca_exn)$sex <- ifelse(colData(sca_exn)$sample_id == "AD2b", "female", 
                               ifelse(colData(sca_exn)$sample_id == "AD4", "male",
                                      ifelse(colData(sca_exn)$sample_id == "C1", "male",
                                             ifelse(colData(sca_exn)$sample_id == "C3", "female", "none!"))))

table(colData(sca_exn)$Group, colData(sca_exn)$sex)


## ----mast_sex_cov2, eval=F, echo=T--------------------------------------------
# zlmCond_exn_sex <- zlm(~ngeneson + sex + Group, sca_exn)
# sumZLM_exn_sex <- summary(zlmCond_exn_sex,
#                   doLRT=paste0("Group","AD"))
# sumDT_exn_sex <- sumZLM_exn_sex$datatable


## ----mast_noRE_sex_saveDT,include=T,echo=F, eval=F----------------------------
# saveRDS(sumDT_exn_sex, "sumDT_exn_sex.rds")


## ----include=T,echo=T, eval=F-------------------------------------------------
# sumDT_exn_sex <- readRDS("~/Downloads/sumDT_exn_sex.rds")


## ----mast_noRE_sex_loadDT,include=T,echo=F, eval=T----------------------------
sumDT_exn_sex <- readRDS("./data/sumDT_exn_sex.rds")


## ----mast_noRE_getRes_sex,include=T,echo=T, eval=T----------------------------
de_mast_exn_sex <- merge(sumDT_exn_sex[sumDT_exn_sex$component=="H"&sumDT_exn_sex$contrast==paste0("Group","AD"),],
               sumDT_exn_sex[sumDT_exn_sex$component=="logFC"&sumDT_exn_sex$contrast==paste0("Group","AD"),],by='primerid')
de_mast_exn_sex <- dplyr::select(de_mast_exn_sex,primerid,coef.y,z.y,`Pr(>Chisq).x`)
names(de_mast_exn_sex) <- c("geneID","log2Fc","z","pvalue")
de_mast_exn_sex$FDR <- p.adjust(de_mast_exn_sex$pvalue,'fdr')
de_mast_exn_sex <- de_mast_exn_sex[order(de_mast_exn_sex$FDR),]

head(de_mast_exn_sex, 3)


## -----------------------------------------------------------------------------
de_mast_exn_sex$sig <- de_mast_exn_sex$FDR < 0.05

ggplot(de_mast_exn_sex, aes(x = log2Fc, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw()


## ----mast_noRE_xist_table,include=T,echo=T, eval=T----------------------------

de_mast_exn_sex[de_mast_exn_sex$geneID == "XIST", ]


## ----mast_noRE_xist_table2,include=T,echo=T, eval=T---------------------------
de_mast_exn[de_mast_exn$geneID == "XIST", ]


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Pseudoreplication bias

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Pseudoreplication bias

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Pseudobulk

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Pseudobulk

---
"    
  )
  
}



## ----pseudo_agg,include=T,echo=T----------------------------------------------
seu_exn_aggr <- Seurat::AggregateExpression(seu_exn, return.seurat = T, group.by = c("sample_id", "condition"))
# get the raw un-normalized counts
counts_aggr <- as.matrix(GetAssayData(seu_exn_aggr, assay = "RNA", layer = "counts"))
head(counts_aggr, 3)


## ----pseudo_prep,include=T,echo=T---------------------------------------------
library(DESeq2)
column_meta <- data.frame(
  row.names = colnames(counts_aggr),
  condition = ifelse(grepl("AD", colnames(counts_aggr)), "AD", "CON")
)
column_meta$condition <- factor(column_meta$condition, levels = c("CON", "AD"))


dds_pseudo_exn <- DESeqDataSetFromMatrix(countData = counts_aggr, colData = column_meta, design = ~condition)
colData(dds_pseudo_exn)


## ----pseudo_filter,include=T,echo=T-------------------------------------------
smallestGroupSize <- 2
keep <- rowSums(counts(dds_pseudo_exn) >= 10) >= smallestGroupSize
table(keep)
dds_pseudo_exn <- dds_pseudo_exn[keep,]


## ----pseudo_results,include=T,echo=T------------------------------------------

dds_pseudo_exn <- DESeq(dds_pseudo_exn)
resultsNames(dds_pseudo_exn)




## ----pseudo_results2,include=T,echo=T-----------------------------------------
res_pseudo_exn <- results(dds_pseudo_exn, name = "condition_AD_vs_CON")
de_pseudo_exn <- res_pseudo_exn %>%
  data.frame %>%
  arrange(pvalue) %>%
  rownames_to_column(var = "geneID") 

head(de_pseudo_exn, 3)


## ----pseudo_MA_PCA, fig.width = 10--------------------------------------------
ma <- ggplot(de_pseudo_exn, aes(x = baseMean, y = log2FoldChange)) + geom_point() + theme_bw() +scale_x_log10() + ggtitle("MA plot") 
pca <- plotPCA(rlog(dds_pseudo_exn))  + coord_cartesian() + theme_bw() + ggtitle("PCA")

ggarrange(ma, pca, ncol = 2)


## ----pseudo_plot,include=T,echo=T, fig.height = 3-----------------------------

de_pseudo_exn$sig = de_pseudo_exn$padj < 0.05
ggplot(de_pseudo_exn %>% na.omit, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw()




## ----pseudo_violin1, echo=T, eval=T, fig.width=7------------------------------
de_pseudo_exn_df <- na.omit(de_pseudo_exn)
top_down_AD_ps <- head(de_pseudo_exn_df[de_pseudo_exn_df$log2FoldChange < 0, ], 5)
VlnPlot(seu_exn, features = top_down_AD_ps$geneID, stack = T, flip = T, pt.size = 1)  + scale_x_discrete(labels=c('CON', 'AD'))



## ----deseq2_sex_cov1, eval=T, echo=T, fig.height = 3--------------------------

sample_ids <- rownames(colData(dds_pseudo_exn))
colData(dds_pseudo_exn)$sex <- ifelse(sample_ids == "AD2b_AD", "female", 
                                      ifelse(sample_ids == "AD4_AD", "male",
                                             ifelse(sample_ids == "C1_CON", "male",
                                                    ifelse(sample_ids == "C3_CON", "female", "none!"))))
colData(dds_pseudo_exn)$sex <- as.factor(colData(dds_pseudo_exn)$sex)
plotPCA(rlog(dds_pseudo_exn), intgroup = "sex")  + coord_cartesian() + theme_bw() + ggtitle("PCA")


## ----deseq2_sex_cov2, eval=T, echo=T------------------------------------------
dds_pseudo_exn_sex <- dds_pseudo_exn
design(dds_pseudo_exn_sex) <- ~sex + condition
dds_pseudo_exn_sex <- DESeq(dds_pseudo_exn_sex)


## ----deseq2_sex_cov3, eval=T, echo=T------------------------------------------

res_pseudo_exn_sex <- results(dds_pseudo_exn_sex, name = "condition_AD_vs_CON")

head(res_pseudo_exn_sex, 3)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Mixed model with random effect

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Mixed model with random effect

---
"    
  )
  
}



## ----mast_re_prep,include=T,echo=T, eval=F------------------------------------
# # recommended here for RE: https://github.com/RGLab/MAST/issues/133
# sca_exn_filt<-filterLowExpressedGenes(sca_exn,threshold=0.1)
# 
# cdr_filt <- colSums(assay(sca_exn_filt)>0)
# colData(sca_exn_filt)$ngeneson <- scale(cdr_filt)
# 


## ----mast_re,include=T,echo=T, eval=F-----------------------------------------
# zlm_re <- zlm(~ngeneson + sex+ Group + (1 | sample_id),
#                sca_exn_filt,
#                ebayes = FALSE,
#                method="glmer")
# 
# sumZLM_re <- summary(zlm_re,
#                   doLRT=paste0("Group","AD"))
# sumDT_re <- sumZLM_re$datatable
# 


## ----mast_re_saveDT,include=T,echo=F, eval=F----------------------------------
# saveRDS(sumDT_re, "sumDT_exn_re.rds")


## ----include=T,echo=T, eval=F-------------------------------------------------
# sumDT_re <- readRDS("~/Downloads/sumDT_exn_re.rds")


## ----mast_re_loadDT,include=T,echo=F, eval=T----------------------------------
sumDT_re <- readRDS("./data/sumDT_exn_re.rds")


## ----mast_re2,include=T,echo=T, eval=T----------------------------------------
de_mast_re <- merge(sumDT_re[sumDT_re$component=="H"&sumDT_re$contrast==paste0("Group","AD"),],
               sumDT_re[sumDT_re$component=="logFC"&sumDT_re$contrast==paste0("Group","AD"),],by='primerid')
de_mast_re <- dplyr::select(de_mast_re,primerid,coef.y,z.y,`Pr(>Chisq).x`)
names(de_mast_re) <- c("geneID","log2Fc","z","pvalue")
de_mast_re$FDR <- p.adjust(de_mast_re$pvalue,'fdr')
de_mast_re <- de_mast_re[order(de_mast_re$FDR),]

head(de_mast_re, 3)


## ----mast_re_plot,include=T,echo=T, eval=T------------------------------------

de_mast_re$sig <- de_mast_re$FDR < 0.05

ggplot(de_mast_re, aes(x = log2Fc, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# DESeq2 and MAST using Seurat FindMarkers?

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# DESeq2 and MAST using Seurat FindMarkers?

---
"    
  )
  
}



## -----------------------------------------------------------------------------
# test deseq output to deseq using findmarkers
Idents(seu_exn_aggr) <- "condition"
seu_deseq <- FindMarkers(object = seu_exn_aggr,
                         ident.1 = "AD",
                         ident.2 = "CON",
                         test.use = "DESeq2", 
                         slot = "counts",
                         min.cells.group = 2)

seu_deseq$sig = seu_deseq$p_val_adj < 0.05


## -----------------------------------------------------------------------------
p_seuD <- ggplot(seu_deseq, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle("FindMarkers")
p_pseudo <- ggplot(de_pseudo_exn %>% na.omit, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle("Using DESeq2 directly")

ggarrange(p_seuD, p_pseudo, common.legend = TRUE, legend="bottom")


## ----eval=F-------------------------------------------------------------------
# cdr_exn <- colSums(seu_exn[["RNA"]]$data>0)
# seu_exn$ngeneson <- scale(cdr_exn)
# 
# Idents(seu_exn_aggr) <- "condition"
# seu_mast <- FindMarkers(object = seu_exn,
#                         group.by = "condition",
#                         ident.1 = "AD",
#                         ident.2 = "CON",
#                         test.use = "MAST",
#                         slot = "data",
#                         latent.vars = "ngeneson"
# )
# 
# seu_mast$sig = seu_mast$p_val_adj < 0.05
# 


## ----saveSeuMast,include=T,echo=F, eval=F-------------------------------------
# saveRDS(seu_mast, "seu_mast_exn.rds")


## ----include=T,echo=T, eval=F-------------------------------------------------
# seu_mast <- readRDS("~/Downloads/seu_mast_exn.rds")


## ----loadSeuMast,include=T,echo=F, eval=T-------------------------------------
seu_mast <- readRDS("./data/seu_mast_exn.rds")


## -----------------------------------------------------------------------------

p_seuM <- ggplot(seu_mast, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle("FindMarkers")
p_mast <- ggplot(de_mast_exn, aes(x = log2Fc, y = -log10(pvalue), color = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle("Using MAST directly")

ggarrange(p_seuM, p_mast, common.legend = TRUE, legend="bottom")


## ----echo=F, warning=FALSE, message=F-----------------------------------------
rm(seu_obj, counts, counts_mat, sca_exn, seu_exn, seu_exn_data,sumDT_exn, sumDT_exn_sex, dds_pseudo_exn,res_pseudo_exn, seu_mast,seu_deseq)
gc()


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# CITE-seq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# CITE-seq

---
"    
  )
  
}



## ----sec3_CITE_prep,include=TRUE,eval=T---------------------------------------

rna.mat <- readRDS("data/pbmc_umi_mtx.rds")
dim(rna.mat)

hto.mat <- readRDS("data/pbmc_hto_mtx.rds")
dim(hto.mat)
hto.mat[1:5,1:5]


## -----------------------------------------------------------------------------
seu_obj <- CreateSeuratObject(counts=rna.mat,project="citeSeq_demo")
seu_obj[["HTO"]] <- CreateAssayObject(counts=hto.mat)
seu_obj


## ----sec3_CITE_clust,include=TRUE,eval=T--------------------------------------

DefaultAssay(seu_obj) <- "RNA"
seu_obj <- data_proc(seu_obj)
seu_obj <- ScaleData(seu_obj)


## -----------------------------------------------------------------------------
seu_obj <- quick_clust(seu_obj)
DimPlot(seu_obj,group.by = "seurat_clusters",pt.size = 0.2,label = TRUE)+NoLegend()


## -----------------------------------------------------------------------------

hto.mat[,1]



## ----echo=F, eval=T, warning=F, message=FALSE, include=F----------------------

rm(hto.mat, rna.mat )
gc()


## ----sec3_CITE_hto,include=TRUE,eval=T----------------------------------------
DefaultAssay(seu_obj) <- "HTO"
seu_obj <- NormalizeData(seu_obj, assay="HTO", normalization.method="CLR")



## -----------------------------------------------------------------------------
seu_obj <- HTODemux(seu_obj, assay = "HTO", positive.quantile = 0.99)

head(seu_obj,2)


## ----CITESeq_posEG,include=TRUE,eval=TRUE,dpi=300-----------------------------
# Distribution of HTO-A level
RidgePlot(seu_obj, features = "HTO-A", group.by = "orig.ident")+NoLegend()



## ----CITESeq_posEG2,include=TRUE,eval=TRUE,dpi=300----------------------------
RidgePlot(seu_obj,
          features = c("HTO-A","HTO-B"),
          group.by = "hash.ID")+NoLegend()


## -----------------------------------------------------------------------------
table(seu_obj$HTO_classification.global)
#
table(seu_obj$hash.ID)
#
table(seu_obj$HTO_classification.global,seu_obj$hash.ID)


## ----CITESeq_spltUMAP,include=TRUE,eval=TRUE,dpi=300--------------------------
DimPlot(seu_obj,group.by = "seurat_clusters",label = TRUE,pt.size = 0.2,split.by = "hash.ID",ncol = 5)+NoLegend()


## -----------------------------------------------------------------------------

seu_obj_onlysinglet <- subset(seu_obj, HTO_classification.global=="Singlet")


seu_obj_onlyHTOA <- subset(seu_obj, HTO_classification=="HTO-A")



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Review

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Review

---
"    
  )
  
}


