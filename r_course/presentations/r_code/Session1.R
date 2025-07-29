params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Introduction to single cell sequencing

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Introduction to single cell sequencing

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Ranger

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger

---
"    
  )
  
}



## wget -O cellranger-8.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.0.tar.gz?Expires=1711772964&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=muvzcbqxba6d-blyYS02MVfLlzwZk6iZNQWXdaoCLnl7owW2nEN-IHwSPwdNoYl-6Xia7rr0S1sLCUQTsekGm2pQKcd0kqK~ndHK0DM7SwSVpXLlRvBV5pXt~EIlsxATVBKVeQLnUy698N-WnRlT~ahjlU-nMdpomX9-lOkF~w8gbgHBdtPXunTWfW87sSJLpHMDVENSF7TFJsXERDwDnsXyQLCuEhfGTCOnupkaATlLEr9kaeCStePKkwGyqgi1m8Ua02NNGHWPIJ6I1mDt695wo~dgptpJF4SDNRTyE-TuXrHfIqRjZB60zhWRJczFo2kpL7FCKwliE-vJ6djcSw__"

## wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

## tar -xzvf cellranger-8.0.0.tar.gz
## tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

## export PATH=/PATH_TO_CELLRANGER_DIRECTORY/cellranger-8.0.0:$PATH

## cellranger count --id=my_run_name \
##    --fastqs=PATH_TO_FASTQ_DIRECTORY \
##    --transcriptome=/PATH_TO_CELLRANGER_DIRECTORY/refdata-gex-GRCh38-2020-A
##    --create-bam=true

## cellranger mkgtf Homo_sapiens.GRCh38.ensembl.gtf \
## Homo_sapiens.GRCh38.ensembl.filtered.gtf \
##                    --attribute=gene_biotype:protein_coding \
##                    --attribute=gene_biotype:lncRNA \
##                    --attribute=gene_biotype:antisense \
##                    --attribute=gene_biotype:IG_LV_gene \
##                    --attribute=gene_biotype:IG_V_gene \
##                    --attribute=gene_biotype:IG_V_pseudogene \
##                    --attribute=gene_biotype:IG_D_gene \
##                    --attribute=gene_biotype:IG_J_gene \
##                    --attribute=gene_biotype:IG_J_pseudogene \
##                    --attribute=gene_biotype:IG_C_gene \
##                    --attribute=gene_biotype:IG_C_pseudogene \
##                    --attribute=gene_biotype:TR_V_gene \
##                    --attribute=gene_biotype:TR_V_pseudogene \
##                    --attribute=gene_biotype:TR_D_gene \
##                    --attribute=gene_biotype:TR_J_gene \
##                    --attribute=gene_biotype:TR_J_pseudogene \
##                    --attribute=gene_biotype:TR_C_gene

## cellranger mkref --genome=custom_reference \
## --fasta=custom_reference.fa  \
## --genes=custom_reference_filtered.gtf

## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Output files

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Output files

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Web Summary QC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Web Summary QC

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Comparison of other scRNA-seq assays

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Comparison of other scRNA-seq assays

---
"    
  )
  
}


