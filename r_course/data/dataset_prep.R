
library(Seurat)

# AD2b <- Read10X_h5("~/Desktop/demo_data/matrix_forcmp/AD2b_cellbender_v0.2.2_filtered.h5")
# AD4 <- Read10X_h5("~/Desktop/demo_data/matrix_forcmp/AD4_cellbender_v0.2.2_filtered.h5")
# C1 <- Read10X_h5("~/Desktop/demo_data/matrix_forcmp/C1_cellbender_v0.2.2_filtered.h5")
# C3 <- Read10X_h5("~/Desktop/demo_data/matrix_forcmp/C3_cellbender_v0.2.2_filtered.h5")

test <- readRDS("~/Downloads/GSE287652_General_clusters.RDS")

AD2b <- subset(test, sample=="AD2b")
AD4 <- subset(test, sample=="AD4")
C1 <- subset(test, sample=="C1")
C3 <- subset(test, sample=="C3")

my_seu <- c(C1, C3, AD2b, AD4)

my_names <- c("Excitatory Neurons",
  "Inhibitory Neurons",
  "OPCs",
  "Pericytes",
  "Oligodendrocytes",
  "Endothelial",
  "Astrocytes",
  "Microglia")

clust_ids <- list(c(0,3,6,8,13,14,16,17,18,19,24,28,30,31,33),
     c(12, 15,25,27,29),
     c(7),
     c(23),
     c(1,2,4,10,11),
     c(26),
     c(5,20,21,32),
     c(9,22))

names(clust_ids) <- my_names

my_seu_filt <- lapply(1:4, function(x){
  seu <- CreateSeuratObject(GetAssayData(my_seu[[x]], assay = "RNA", layer="counts"), min.features =200, min.cells=10 )
  seu$percent.mt <- my_seu[[x]]$percent.mt
  seu$scrubletScore <- my_seu[[x]]$scrubletScore
  seu$sample_id <- my_seu[[x]]$sample
  seu$paper_cluster <- my_seu[[x]]$seurat_clusters
  seu <- NormalizeData(seu, normalization.method="LogNormalize")
  #seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- FindVariableFeatures(seu, select.method="vst", nfeatures=3000)
  seu <- ScaleData(seu)
  #seu_filt <- subset(seu, percent.mt < 2.5)
  seu_filt <- ScaleData(seu, vars.to.regress = "percent.mt")
  DefaultAssay(seu_filt) <- "RNA"
  seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)
  seu_filt <- FindNeighbors(seu_filt, reduction = "pca")
  seu_filt <- RunUMAP(seu_filt, dims=1:30, reduction = "pca")
  #seu_filt$sample_id <-  c("C1", "C3", "AD2b", "AD4")[x]
  seu_filt$group <-  c("C", "C", "AD", "AD")[x]
  
  # ids <- test$seurat_clusters[match(paste0(Cells(seu_filt),"_",c("AD2b", "AD4", "C1", "C3")[x]), names(test$seurat_clusters))]
  # names(ids) <- Cells(seu_filt)
  # seu_filt$paper_cluster <- ids
  
  seu_filt$paper_annot <- NA
  
  for(z in 1:length(clust_ids)){
    seu_filt$paper_annot[seu_filt$paper_cluster %in% clust_ids[[z]]] <-  names(clust_ids)[z]
  }
  
  return(seu_filt)
})


my_seu_list <- my_seu_filt

saveRDS(my_seu_filt,"~/Desktop/to_integrate.rds" )
saveRDS(my_seu_filt, "scRNASeq/inst/extdata/data/to_integrate.rds")
#my_seu_list <- readRDS("scRNASeq/inst/extdata/data/to_integrate.rds")


library(Matrix)
mat <- read.csv("~/Downloads/matrix.csv", row.names = 1)
mat2 <- as.matrix(t(mat))
mat3 <- Matrix(mat2,sparse=TRUE)
mat2 <- NULL
mat <- NULL
meta <- read.csv("~/Downloads/metadata (1).csv")
seu <- CreateSeuratObject(mat3, min.features =200, min.cells=10 )
seu <- NormalizeData(seu, normalization.method="LogNormalize")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- FindVariableFeatures(seu, select.method="vst", nfeatures=3000)
seu <- ScaleData(seu)
seu_filt <- seu
#seu_filt <- subset(seu, percent.mt < 2.5)
#seu_filt <- ScaleData(seu, vars.to.regress = "percent.mt")
DefaultAssay(seu_filt) <- "RNA"
seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)
seu_filt <- FindNeighbors(seu_filt, reduction = "pca")
seu_filt <- RunUMAP(seu_filt, dims=1:30, reduction = "pca")

seu_filt$class_label <- meta$class_label
seu_filt$cluster_label <- meta$cluster_label
seu_filt$subclass_label <- meta$subclass_label
seu_filt$cell_type_alias_label <- meta$cell_type_alias_label

saveRDS(seu_filt,"~/Desktop/abm.rds" )
saveRDS(seu_filt, "scRNASeq/inst/extdata/data/abm.rds")


files <- c("~/Downloads/FCAImmP7277552/GRCh38",
  "~/Downloads/FCAImmP7277553/GRCh38",
  "~/Downloads/FCAImmP7277560/GRCh38",
  "~/Downloads/FCAImmP7277561/GRCh38")

labels <- c("liver_CD45pos",
            "liver_CD45neg",
            "liver_CD45pos",
            "liver_CD45neg")

seu_list <- lapply(1:4, function(x){

mtx <- Read10X(files[x])
seu <- CreateSeuratObject(mtx, min.features =200, min.cells=10 )
seu[["dset"]] <- labels[x]# Create a category for sample
seu <- Seurat::RenameCells(seu, add.cell.id=labels[x]) # add sample name in front of cell 
seu <- NormalizeData(seu, normalization.method="LogNormalize")
seu <- FindVariableFeatures(seu, select.method="vst", nfeatures=3000)
seu <- ScaleData(seu)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
feat_s <- cc.genes$s.genes
feat_g2m <- cc.genes$g2m.genes
seu<- CellCycleScoring(seu, s.features = feat_s, g2m.features = feat_g2m)
seu_filt <- subset(seu, percent.mt < 10)
seu_filt <- ScaleData(seu, vars.to.regress = c("percent.mt","Phase"))
DefaultAssay(seu_filt) <- "RNA"
seu_filt <- RunPCA(seu_filt, assay = "RNA", npcs = 50)
seu_filt <- FindNeighbors(seu_filt, reduction = "pca")
seu_filt <- RunUMAP(seu_filt, dims=1:30, reduction = "pca")

return(seu_filt)
})
