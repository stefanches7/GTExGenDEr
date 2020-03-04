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
library(RColorBrewer)
library(tweeDEseqCountData)
library(gplots)

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Overview/globalOverview.R")


dt_tpm <- fread(file = snakemake@input[["countsTpm"]], fill = T, skip = 2)
counts <- fread(file = snakemake@input[["counts"]], fill = T, skip = 2)

#' ## Annotation

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% colnames(dt_tpm)]
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

#samples by sex plot

tissueBySex <- anno[, .N, by = list(SMTSD, SEX)]
decreasingSampleNo <- anno[,.N, by = SMTSD][order(N, decreasing = T),]
tissueBySex[,SMTSD:=factor(SMTSD, levels = unique(decreasingSampleNo$SMTSD))]
ggplot(data = decreasingSampleNo) +
  geom_text(aes(x=factor(SMTSD, levels = SMTSD),y=N,label=rep("", length(N)-4),tail(as.character(N),4)),vjust=0) +
  geom_bar(data = tissueBySex, aes(x = SMTSD, y = N, fill = SEX), stat = "identity") +
  labs(title = "Different tissues in genders",x = "", y = "Count") + 
  theme_light() +
  theme(
    plot.margin = margin(0, 1, 1.2, 1, "cm"), 
    axis.text.x = element_text(face="plain", 
                               size=10, angle=90, hjust = 1, vjust = 0,
                               margin = margin(0, 0, 0, 0, "cm")))


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

deInTissues <- dt_dt[, .N, by = "Tissue"][order(N, decreasing = T)]
ggplot(deInTissues, aes(x = factor(Tissue, levels = unique(deInTissues$Tissue)), y = N)) + 
  geom_bar(stat = "identity",fill ="orange1") +
  scale_y_continuous(breaks = round(seq(0, 14000, by = 2000),1)) +
  labs(title = "Differential genes in tissues",x = "", y = "Count") + 
  theme_light() +
  theme(
    plot.margin = margin(0, 1, 1.2, 1, "cm"), 
    axis.text.x = element_text(face="plain", 
                               size=10, angle=90, hjust = 1, vjust = 0,
                               margin = margin(0, 0, 0, 0, "cm")))

tissuesList <- unique(dt_dt$Tissue)
tissuesDEGenes <- lapply(tissuesList, function(tissue) dt_dt[Tissue == tissue, id])
names(tissuesDEGenes) <- tissuesList
similarityMtx <- matrix(nrow = length(tissuesList), ncol = length(tissuesList))
tissueSimilarityMetric <- function(tissueA, tissueB) {
  length(intersect(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]])) / 
    length(union(tissuesDEGenes[[tissueA]], tissuesDEGenes[[tissueB]]))
}
for (i in 1:length(tissuesList)){
  for (j in 1:length(tissuesList)){
    if (i == j) {
      similarityMtx[i,j] <- 0
      next
      }
    similarityMtx[i,j] <- tissueSimilarityMetric(tissuesList[[i]], tissuesList[[j]])
  }
}
colnames(similarityMtx) <- tissuesList
rownames(similarityMtx) <- tissuesList
heatmap.2(similarityMtx, trace = "none", col = brewer.pal(11, name = "RdBu"), main = "Correlation of tissue DE genes", 
          key.title = "Spearman ρ", margin = c(12, 10))

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

#' ### Sex-specific genes

data("genderGenes")
geneEnsIds <- tstrsplit(dt_tpm$Name, "\\.")[[1]]
Ymale <- geneEnsIds %in% msYgenes
Xescape <- geneEnsIds %in% XiEgenes
maleSamples <- anno[SEX == "M", SAMPID]
femaleSamples <- anno[SEX == "F", SAMPID]

sexCountP <- data.table(geneName = character(), lpval = numeric(), minavg = numeric(), gene = character())
for (maleIdx in which(Ymale)){
  geneName <- counts$Name[maleIdx]
  for (tissue in tissuesList) {
    if (!(geneName %in% dt_dt[Tissue == tissue, id])) next
    pval <- dt_dt[Tissue == tissue & id == geneName, FDR]
    avgMale <- counts[Name == geneName, ..maleSamples] %>% as.numeric %>% mean
    avgFemale <- counts[Name == geneName, ..femaleSamples] %>% as.numeric %>% mean
    sexCountP <- rbindlist(list(sexCountP, list(geneName, -log10(pval), min(avgMale, avgFemale), "Y male")))
  }
}
for (femaleIdx in which(Xescape)){
  geneName <- counts$Name[femaleIdx]
  for (tissue in tissuesList) {
    if (!(geneName %in% dt_dt[Tissue == tissue, id])) next
    pval <- dt_dt[Tissue == tissue & id == geneName, FDR]
    avgMale <- counts[Name == geneName, ..maleSamples] %>% as.numeric %>% mean
    avgFemale <- counts[Name == geneName, ..femaleSamples] %>% as.numeric %>% mean
    sexCountP <- rbindlist(list(sexCountP, list(geneName, -log10(pval), min(avgMale, avgFemale), "X escaping")))
  }
}
sexCountP[,sex:=factor(sex)]
ggplot(sexCountP, aes(x = lpval, y = minavg, color = sex)) + geom_point() + 
  labs(title = "Gender-specific gene expression",x = "-log10(FDR) in DE analysis", y = "Min average TPM (males/females)") + 
  theme_light()
  

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
