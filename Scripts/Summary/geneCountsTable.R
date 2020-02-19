#'---
#' title: Compare brain cortex
#' author: Stefan Dvoretskii
#' wb:
#'  input:
#'  - counts: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
#'  - sampleAnno: "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#'  - phenotypeAnno: "data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/geneCountsTable.R")


setwd("~/projects/GTExGenDEr/")
library(data.table)
library(ggplot2)
library(MASS)
library(magrittr)
library(edgeR)
library(tweeDEseqCountData)
library(DT)

#' ## Read in the data

#dt <- fread(cmd = paste0("zcat < ", snakemake@input[["counts"]]))
dt_tpm <- fread(cmd = paste0("zcat < ", snakemake@input[["countsTpm"]]))

#' ## Annotation

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% colnames(dt_tpm)]
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

#' How many males(=1) vs. females(=2)?
table(anno$SEX)

#' How many samples by subject?
hist(anno[, .N, by = "SUBJID"]$N, main = "Samples by subject", xlab = "# samples", xlim = c(0, 100), breaks = 100)
anno[, .N, by = "SUBJID"][N>100,]

#' What's the difference between gene counts and gene tpm datasets?
dt_tpm[2,5]
#' -> difference in file size is only because of integer vs. float capacity

#' ## Gender sample differences 
#' How many samples for each tissue in male vs. female? 

tissueBySex <- anno[, .N, by = list(SMTSD, SEX)]
ggplot(tissueBySex, aes(x = SMTSD, y = N, fill = SEX)) + geom_bar(stat = "identity") +
  labs(title = "Different tissues in genders",x = "", y = "Count") + 
  theme_light() +
  theme(
    plot.margin = margin(0, 1, 1.2, 1, "cm"), 
    axis.text.x = element_text(face="plain", 
                                 size=7, angle=45, hjust = 1, vjust = 0.5,
                                 margin = margin(-1.7, 0, 0, 0, "cm")))

#' Which tissues do only females have? 
setdiff(tissueBySex[SEX=="F",SMTSD], tissueBySex[SEX=="M",SMTSD])
setdiff(tissueBySex[SEX=="M",SMTSD], tissueBySex[SEX=="F",SMTSD])

nrow(anno[is.na(SMTSD),]) == 0
commonTissuesAnno <- anno[SMTSD %in% 
                            intersect(tissueBySex[SEX=="F",SMTSD], tissueBySex[SEX=="M",SMTSD]), ]

#' ### Analysis of brain cortex differences

maleTissIdx <- anno[SMTSD == "Brain - Cortex" & SEX == "M", SAMPID]
femTissIdx <- anno[SMTSD == "Brain - Cortex" & SEX == "F", SAMPID]

X <- cbind(dt_tpm[, ..femTissIdx], dt_tpm[, ..maleTissIdx]) %>% as.matrix
# TODO remove all zero rows

#' fit glm with edgeR

conditions <- c(rep(0, length(femTissIdx)), rep(1, length(maleTissIdx)))

data(annotEnsembl63)
geneAnno <- annotEnsembl63[,c("Symbol","Chr")] # TODO change to more recent/biomaRt?

geneEnsIds <- tstrsplit(dt_tpm$Name,"\\.")[[1]]
y <- DGEList(counts = X, group = conditions, genes = geneAnno[geneEnsIds,])
design <- model.matrix(~conditions) # design of DE analysis

y <- estimateDisp(y, design = design)
plotBCV(y)
fit <- glmQLFit(y, design = design)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
tt <- topTags(qlf, adjust.method = "BH", p.value = 0.05, n = nrow(counts))

deGeneRows <- rownames(tt) %>% as.numeric

diffGenes <- dt_tpm[deGeneRows, Name]
diffTable <- data.table(tt$table)[,ENS_ID:=diffGenes]
cols <- colnames(diffTable)[3:5]
diffTable[,(cols):=round(.SD, 5), .SDcols = cols]
#DT::datatable(diffTable) # visual table for the server

summary(decideTests(qlf))
plotMD(qlf)
abline(h = c(-1, 1), col = "blue")

#' get enriched GO/KEGG categories with limma
go <- goana(qlf)
goCat <- topGO(go, ont="BP", sort="Up", n=100, truncate=100) %>% as.data.table
cols <- c("P.Up", "P.Down")
#DT::datatable(goCat[P.Up < 0.05 | P.Down < 0.05,(cols):=round(.SD, 5), .SDcols = cols])
kegg <- kegga(qlf)
keggCat <- topKEGG(kegg, sort = "Up")
data(genderGenes)

#' as it seems no big differential hits on categories

#' ## Compare DE genes to already known Y-specific/X-escaping(=female specific) genes
data("genderGenes")
Ymale <- geneEnsIds %in% msYgenes
Xescape <- geneEnsIds %in% XiEgenes

index <- list(Y=Ymale, X=Xescape)
fry(y, index=index, design = design) # limma gene set rotation tests
barcodeplot(qlf$table$logFC, index[[1]], index[[2]])
#' red is women, blue is men

with(qlf$table, plot(logCPM,logFC,pch=16,cex=0.2))
with(qlf$table, points(logCPM[Ymale],logFC[Ymale],pch=16,col="red"))
with(qlf$table, points(logCPM[Xescape],logFC[Xescape],pch=16,col="dodgerblue"))
legend("bottomleft",legend=c("Ymale genes","Xescape genes"),
                                                         pch=16,col=c("red","dodgerblue"))

