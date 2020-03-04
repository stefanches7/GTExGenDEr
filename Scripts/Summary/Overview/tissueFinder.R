#'---
#' title: Find tissue proximity based on the DE genes
#' author: Stefan Dvoretskii
#' wb:
#'  py:
#'  - | 
#'   with open('data/generalTissues.txt') as f:
#'    gt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'   with open('data/detailedTissues.txt') as f:
#'    dt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'  input:
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
#'  - counts: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
#'  - sampleAnno: "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#'  - phenotypeAnno: "data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#'  - detailedTissuesDE: "`sm expand('{wbPD}/DetailedTissues/{tissue}/deTopHits.tsv', tissue = dt)`"
#'  - detailedTissuesGO: "`sm expand('{wbPD}/DetailedTissues/{tissue}/goTopHits.tsv', tissue = dt)`"
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

knitr::opts_chunk$set(fig.width = 14, fig.height = 10)
library(data.table)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(tweeDEseqCountData)

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Overview/globalOverview.R")


detailedTissuesResults <- paste0(snakemake@wildcards[["wbPD"]], "/DetailedTissues/")
generalTissuesResults <- paste0(snakemake@wildcards[["wbPD"]], "/GeneralTissues/")

dt_dt <- NULL
for (tissueDir in dir(detailedTissuesResults)) {
  tableFile <- paste0(detailedTissuesResults, tissueDir, "/deTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(dt_dt)) { dt_dt <-  tissueDEResults}
  else { dt_dt <- rbind(dt_dt, tissueDEResults)}
}

tissuesList <- setdiff(unique(dt_dt$Tissue), c("Lung"))
tissuesDEGenes <- lapply(tissuesList, function(tissue) dt_dt[Tissue == tissue, id])
names(tissuesDEGenes) <- tissuesList
similarityMtx <- matrix(nrow = length(tissuesList), ncol = length(tissuesList))
tissueSimilarityMetric <- function(tissueA, tissueB) {
  length(intersect(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]])) / 
    length(union(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]]))
}
for (i in 1:length(tissuesList))
  for (j in 1:length(tissuesList))
    similarityMtx[i,j] <- tissueSimilarityMetric(tissuesList[[i]], tissuesList[[j]])

colnames(similarityMtx) <- tissuesList
heatmap.2(similarityMtx, trace = "none", col = brewer.pal(11, name = "RdBu"), main = "Correlation of tissue DE genes of 45 known tissues", 
          key.title = "Spearman ρ", margin = c(12, 10))

tableFile <- paste0(detailedTissuesResults, "Lung", "/deTopHits.tsv")
tissueXDEResults <- fread(tableFile)
#' ### Tissue X DE results 
DT::datatable(tissueXDEResults[,Tissue:=NULL])

#' ### Tissue X results

tissuesList <- unique(dt_dt$Tissue)
tissuesList <- sapply(tissuesList, function(tiss) ifelse(tiss == "Lung", "X", tiss))
tissuesDEGenes <- lapply(tissuesList, function(tissue) dt_dt[Tissue == tissue, id])
names(tissuesDEGenes) <- tissuesList
tissuesDEGenes[["X"]] <- dt_dt[Tissue == "Lung", id]
similarityMtx <- matrix(nrow = length(tissuesList), ncol = length(tissuesList))
tissueSimilarityMetric <- function(tissueA, tissueB) {
  length(intersect(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]])) / 
    length(union(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]]))
}
for (i in 1:length(tissuesList))
  for (j in 1:length(tissuesList))
    similarityMtx[i,j] <- tissueSimilarityMetric(tissuesList[[i]], tissuesList[[j]])

colnames(similarityMtx) <- tissuesList
note <- matrix("", nrow=nrow(c), ncol=ncol(c))
for (i in 1:nrow(c)) {
  note[which(tissuesList == "X"), i] <- "*"  
  note[i, which(tissuesList == "X")] <- "*"  
}

heatmap.2(similarityMtx, trace = "none", col = brewer.pal(11, name = "RdBu"), main = "Correlation of tissue DE genes with tissue X", 
          key.title = "Spearman ρ", margin = c(12, 10), cellnote = note, notecol = "red")

