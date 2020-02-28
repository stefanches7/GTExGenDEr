#'---
#' title: Comparison of males and females in a single tissue
#' author: Stefan Dvoretskii
#' wb:
#'  input:
#'  - scriptMappingFlag: "Output/scriptMapping.done"
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
#'  - sampleAnno: "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#'  - phenotypeAnno: "data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#'  - geneAnno: "data/geneAnno.tsv"
#'  output:
#'  - deTopHits: "{wbPD_PP}/deTopHits.tsv"
#'  - goTopHits: "{wbPD_PP}/goTopHits.tsv"
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
source("Scripts/_Template/utils.R")
parseWBHeader("Scripts/_Template/singleTissueAnalysis.R")

tissueType <- snakemake@wildcards[["wbP"]]

currentTissue <- recoverTissueName(snakemake@wildcards[["wbPP"]])
print(paste0("For tissue ", currentTissue))

library(data.table)
library(ggplot2)
library(magrittr)
library(tweeDEseqCountData)
library(edgeR)
library(DT)

#' ## Read in the data

dt_tpm <- fread(file = snakemake@input[["countsTpm"]], fill = T, skip = 2)

#' ## Annotation

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% colnames(dt_tpm)]
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

#' ## Analysis of differences

if (tissueType == "GeneralTissues") {
  maleTissIdx <- anno[SMTS == currentTissue & SEX == "M", SAMPID]
  femTissIdx <- anno[SMTS == currentTissue & SEX == "F", SAMPID]
} else {
  maleTissIdx <- anno[SMTSD == currentTissue & SEX == "M", SAMPID]
  femTissIdx <- anno[SMTSD == currentTissue & SEX == "F", SAMPID]
}
X <- cbind(dt_tpm[, ..femTissIdx], dt_tpm[, ..maleTissIdx]) %>% as.matrix
# TODO remove all zero rows

#' ### GLM

conditions <- c(rep(0, length(femTissIdx)), rep(1, length(maleTissIdx)))

geneAnno <- fread(snakemake@input[["geneAnno"]])

geneEnsIds <- dt_tpm$Name

y <- DGEList(counts = X, group = conditions, genes = geneAnno)

design <- model.matrix(~conditions) # design of DE analysis

y <- estimateDisp(y, design = design)
plotBCV(y)
fit <- glmQLFit(y, design = design)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)

#' ### Differential genes
tt <- topTags(qlf, adjust.method = "BH", p.value = 0.05, n = nrow(X))

deGeneRows <- rownames(tt) %>% as.numeric

diffGenes <- dt_tpm[deGeneRows, Name]
diffTable <- data.table(tt$table)[,ENS_ID:=diffGenes]

chromosomeCount <- table(diffTable$chromosome_name)
t <- sort(chromosomeCount, decreasing = T)[1:10]
plt <- data.table(name = factor(names(t), levels = names(t)), value = t %>% as.numeric())
ggplot(plt, aes(name, value)) + geom_bar(stat = "identity", fill = "springgreen2") +
  labs(title = "Top-10 differential chromosomes", xlab = "Chromosome", ylab = "Differential genes count") + theme_light()

diffTable[,Tissue:=currentTissue]
write.table(diffTable, file = snakemake@output[["deTopHits"]], sep = "\t", row.names = F)


cols <- c("logFC", "logCPM", "F")
diffTable[,(cols):=round(.SD, 5), .SDcols = cols]
DT::datatable(diffTable) # visual table for the server

summary(decideTests(qlf))
plotMD(qlf)
abline(h = c(-1, 1), col = "blue")

#' ### Enriched GO categories
go <- goana(qlf)
goUp <- topGO(go, ont="BP", sort="Up", number=200) %>% as.data.table
goDown <- topGO(go, ont="BP", sort="Down", number=200) %>% as.data.table

goCat <- rbind(goUp, goDown)

cols <- c("P.Up", "P.Down")
DT::datatable(goCat[P.Up < 0.05 | P.Down < 0.05,(cols):=round(.SD, 5), .SDcols = cols])
goCat[,Tissue:=currentTissue]
write.table(goCat, file = snakemake@output[["goTopHits"]], sep = "\t", row.names = F)
#kegg <- kegga(qlf)
#keggCat <- topKEGG(kegg, sort = "Up")


#' ## Compare DE genes to already known Y-specific/X-escaping(=female specific) genes

data("genderGenes")
geneEnsIds <- tstrsplit(dt_tpm$Name, "\\.")[[1]]
Ymale <- geneEnsIds %in% msYgenes
Xescape <- geneEnsIds %in% XiEgenes

index <- list(Y=Ymale, X=Xescape)

fry(y$counts, index=index, design = design) # limma gene set rotation tests
barcodeplot(qlf$table$logFC, index[[1]], index[[2]])
#' red is women, blue is men

with(qlf$table, plot(logCPM,logFC,pch=16,cex=0.2))
with(qlf$table, points(logCPM[Ymale],logFC[Ymale],pch=16,col="red"))
with(qlf$table, points(logCPM[Xescape],logFC[Xescape],pch=16,col="dodgerblue"))
legend("bottomleft",legend=c("Ymale genes","Xescape genes"),
       pch=16,col=c("red","dodgerblue"))

#' ### Differential genes excluding apart from X, Y chromosomes
DT::datatable(data.table(tt$table)[chromosome_name != "X" & chromosome_name != "Y",])
