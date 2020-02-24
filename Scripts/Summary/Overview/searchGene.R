#'---
#' title: Lookup tables for specific gene/GO term across tissues
#' author: Stefan Dvoretskii
#' wb:
#'  py:
#'  - | 
#'   with open('data/generalTissues.txt') as f:
#'    gt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'   with open('data/detailedTissues.txt') as f:
#'    dt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'  input:
#'  - detailedTissuesDE: "`sm expand('{wbPD}/DetailedTissues/{tissue}/deTopHits.tsv', tissue = dt)`"
#'  - detailedTissuesGO: "`sm expand('{wbPD}/DetailedTissues/{tissue}/goTopHits.tsv', tissue = dt)`"
#'  - generalTissuesDE: "`sm expand('{wbPD}/GeneralTissues/{tissue}/deTopHits.tsv', tissue = gt)`"
#'  - generalTissuesGO: "`sm expand('{wbPD}/GeneralTissues/{tissue}/goTopHits.tsv', tissue = gt)`"
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

library(data.table)
library(ggplot2)

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Overview/searchGene.R")


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
gt_dt <- NULL
for (tissueDir in dir(generalTissuesResults)) {
  tableFile <- paste0(generalTissuesResults, tissueDir, "/deTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(gt_dt)) { gt_dt <-  tissueDEResults}
  else { gt_dt <- rbind(gt_dt, tissueDEResults)}
}

#' ## All DE genes in all tissues
DT::datatable(rbind(dt_dt, gt_dt) %>% unique)

go_dt <- NULL
for (tissueDir in dir(detailedTissuesResults)) {
  tableFile <- paste0(detailedTissuesResults, tissueDir, "/goTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(go_dt)) { go_dt <-  tissueDEResults}
  else { go_dt <- rbind(go_dt, tissueDEResults)}
}
names(go_dt)[3] <- "Size"

go_gt <- NULL
for (tissueDir in dir(generalTissuesResults)) {
  tableFile <- paste0(generalTissuesResults, tissueDir, "/goTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(go_gt)) { go_gt <-  tissueDEResults}
  else { go_gt <- rbind(go_gt, tissueDEResults)}
}
names(go_gt)[3] <- "Size"

DT::datatable(rbind(go_dt, go_gt) %>% unique)

