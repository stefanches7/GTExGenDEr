#'---
#' title: Load genetic annotation
#' author: Stefan Dvoretskii
#' wb:
#'  input:
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
#'  output:
#'  - geneAnno: "data/geneAnno.tsv"
#'  type: script
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

library(data.table)
library(biomaRt)

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Structure/loadGeneAnno.R")

dt_tpm <- fread(file = snakemake@input[["countsTpm"]], fill = T, skip = 2)

ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
anno <- getBM(attributes = c("description", "chromosome_name", "hgnc_symbol", "ensembl_gene_id_version"), 
      filters = "ensembl_gene_id_version", values = dt_tpm$Name, mart = ensembl)
a <- anno[,1:3]
rownames(a) <- anno$ensembl_gene_id_version
write.table(anno, file = snakemake@output[["geneAnno"]], row.names = F, sep = "\t")