#'---
#' title: Differential expression/GO enrichment overview across tissues
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
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---


library(data.table)
library(ggplot2)

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Overview/globalOverview.R")

dt_dt <- NULL
detailedTissuesResults <- paste0(snakemake@wildcards[["wbPD"]], "/DetailedTissues/")
generalTissuesResults <- paste0(snakemake@wildcards[["wbPD"]], "/GeneralTissues/")

for (tissueDir in dir(detailedTissuesResults)) {
  tableFile <- paste0(detailedTissuesResults, tissueDir, "/deTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(dt_dt)) { dt_dt <-  tissueDEResults}
  else { dt_dt <- rbind(dt_dt, tissueDEResults)}
}
go_dt <- NULL

for (tissueDir in dir(detailedTissuesResults)) {
  tableFile <- paste0(detailedTissuesResults, tissueDir, "/goTopHits.tsv")
  if (!file.exists(tableFile)) next
  tissueDEResults <- fread(tableFile)
  if (is.null(go_dt)) { go_dt <-  tissueDEResults}
  else { go_dt <- rbind(go_dt, tissueDEResults)}
}
names(go_dt)[3] <- "Size"

symbolsCount <- table(dt_dt$hgnc_symbol)
sort(symbolsCount, decreasing = T)[1:10]

ensgCount <- table(dt_dt$id)
sort(ensgCount, decreasing = T)[1:10]

#' ## Overall DE analysis

hist(symbolsCount[2:length(symbolsCount)] %>% as.numeric, breaks = 50, main = "Gene symbols in tissues", 
     xlab = "# occured across tissues", ylab = "Count gene symbols")
hist(ensgCount[2:length(symbolsCount)] %>% as.numeric, breaks = 50, main = "Genes in tissues", 
     xlab = "# occured across tissues", ylab = "Count gene id")

sum(symbolsCount > 15)
hist(symbolsCount[2:length(symbolsCount)] %>% as.numeric, breaks = 50, xlim = c(0,10), xlab = "Times occured significant", 
     main = "Significance of genes count across general tissues")

#hist(symbolsCountGen %>% as.numeric, breaks = 10, main = "Significance of genes count across general tissues",
#     xlab = "Times occured significant")

chromosomeCount <- table(dt_dt$chromosome_name)
t <- sort(chromosomeCount, decreasing = T)[1:10]
plt <- data.table(name = factor(names(t), levels = names(t)), value = t %>% as.numeric())
ggplot(plt, aes(name, value)) + geom_bar(stat = "identity", fill = "springgreen2") +
  labs(title = "Top-10 differential chromosomes", xlab = "Chromosome", ylab = "Differential genes count") + theme_light()

#' ## Overall GO analysis
termCount <- table(go_dt$Term)
sort(termCount, decreasing = T)[1:10]
#' ### GO terms times enriched in detailed tissues
dt <- merge(go_dt[, .(Term, Size)], go_dt[, .N, by = "Term"])
colnames(dt) <- c("Term", "Size", "# occured across tissues")
ggplot(dt, aes(Size, `# occured across tissues`)) + geom_point()
tissues <- vapply(dt$Term, function(t) paste(go_dt[Term == t,Tissue], collapse = ",", sep = ",\n"), FUN.VALUE = character(1)) %>%
  as.character
dt[,Tissues:=tissues]
DT::datatable(unique(dt))
