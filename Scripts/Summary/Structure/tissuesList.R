#'---
#' title: Get available tissues from annotation
#' author: Stefan Dvoretskii
#' wb:
#'  input:
#'  - counts: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
#'  - sampleAnno: "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#'  - phenotypeAnno: "data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#'  output:
#'  - generalTissues: "data/generalTissues.txt"
#'  - detailedTissues: "data/detailedTissues.txt"
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Tissues/Single_Tissues/tissuesList.R")

library(data.table)
library(magrittr)

setwd("~/projects/GTExAna/")

con <- gzfile(snakemake@input[["countsTpm"]])

header <- readLines(con, n = 3)[[3]]
availableSamples <-  strsplit(header, "\t")[[1]]

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% availableSamples,]
unique(samples[SMTS == "Blood", SMTSD])
unique(samples[SMTS == "Blood Vessel", SMTSD])
unique(samples[SMTS == "Brain", SMTSD])
#' 
#' no real **aggregation** needed - each sample already contains the overall tissue

samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

commonSMTS <- intersect(anno[SEX=="F",SMTS], anno[SEX=="M",SMTS]) # "general" (e.g. brain)
commonSMTSD <- intersect(anno[SEX=="F",SMTSD], anno[SEX=="M",SMTSD]) # "detailed" (e.g. Brain - Cortex) 

fh <- file(snakemake@output[["generalTissues"]])
writeLines(commonSMTS, fh)
close(fh)
fh <- file(snakemake@output[["detailedTissues"]])
writeLines(commonSMTSD, fh)
close(fh)