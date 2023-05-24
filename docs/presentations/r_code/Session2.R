params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(bioMart))
suppressPackageStartupMessages(require(SeuratWrappers))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----load_pack,include=TRUE,eval=FALSE----------------------------------------
## library(Seurat)
## library(ggplot2)
## library(scran)
## library(scuttle)
## library(reticulate)


## -----------------------------------------------------------------------------
download.file("https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz","10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")

# untar("10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.tar.gz")



