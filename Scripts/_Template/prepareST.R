#'---
#' title: Prepare single-tissue analyses 
#' author: Stefan Dvoretskii
#' wb:
#'  py:
#'  - | 
#'   with open('data/generalTissues.txt') as f:
#'    gt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'   with open('data/detailedTissues.txt') as f:
#'    dt = [line.replace("(", "8_").replace(")", "_9").replace("-","_").replace(" ", "_") for line in f.read().splitlines()]
#'  input:
#'  - generalTissues: "data/generalTissues.txt"
#'  - detailedTissues: "data/detailedTissues.txt"
#'  - countsTpm: "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
#'  - sampleAnno: "data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#'  - phenotypeAnno: "data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#'  output:
#'  - countsGen: "`sm expand(Output/ProcessedData/GeneralTissues/{tissue}/counts.txt, tissue=gt)"
#'  - countsDet: "`sm expand(Output/ProcessedData/DetailedTissues/{tissue}/counts.txt, tissue=dt)"
#'  type: noindex
#'---

setwd("~/projects/GTExGenDEr/")
source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/Summary/Structure/prepareST.R")

dt_tpm <- fread(file = snakemake@input[["countsTpm"]], fill = T, skip = 2)

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% colnames(dt_tpm)]
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

#' ### Analysis of differences
if (tissueType == "GeneralTissues") {
  maleTissIdx <- anno[SMTS == currentTissue & SEX == "M", SAMPID]
  femTissIdx <- anno[SMTS == currentTissue & SEX == "F", SAMPID]
} else {
  maleTissIdx <- anno[SMTSD == currentTissue & SEX == "M", SAMPID]
  femTissIdx <- anno[SMTSD == currentTissue & SEX == "F", SAMPID]
}
X <- cbind(dt_tpm[, ..femTissIdx], dt_tpm[, ..maleTissIdx]) %>% as.matrix
