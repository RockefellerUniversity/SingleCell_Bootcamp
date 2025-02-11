params <-
list(isSlides = "yes")

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



## wget -O cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1701688001&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=Puwpsqsf~wMQz5e~PwTvM2DRQO1XdJ~9zeLCWqX6tVbOx~dnf24hP1bwlmNhybr3SZUQ8C12ywcICMH6Au02wxiCRm1uuTxZ0Uvq8g~s8L8s6XFyhepdi6Qjq8dzXNGoxswg3hModjKWVptTWq-MTHBDZv~yTFB7QAM9lzHHXo6SPWg8Fnx30ngmtGC5tDReVOiJ3DY0hsFvZtG3HaQ-HEbnzEH3lre-0rpWMBlsQu-vZ4RnKE0o3Xv6pQsb6261M19nHcpCsGhDCkFjDDbradx~SNw5rpY-HMxM4SnRuaOOI0rYyDNn7xdTat3eFj7rlgATXRaYx7SYNqDYKSrNWw__"

## wget -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

## tar -xzvf cellranger-7.2.0.tar.gz
## tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

## export PATH=/PATH_TO_CELLRANGER_DIRECTORY/cellranger-7.1.0:$PATH

## cellranger count --id=my_run_name \
##    --fastqs=PATH_TO_FASTQ_DIRECTORY \
##    --transcriptome=/PATH_TO_CELLRANGER_DIRECTORY/refdata-gex-GRCh38-2020-A

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

# Cell Ranger -  Output files

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger - Output files

---
"    
  )
  
}



## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Cell Ranger - Web Summary QC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Cell Ranger - Web Summary QC

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


