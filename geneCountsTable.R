setwd("~/projects/GTExAna/")
library(data.table)

#' ## Read in the data

dt <- fread("zcat < data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
dt_tpm <- fread("zcat < data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

#' ## Annotation

samples <- fread("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] #take sample id, read only characters before the *second* horizontal dash
phenotypes <- fread("data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)

#' How many males(=1) vs. females(=2)?
table(anno$SEX)

#' How many samples by subject?
hist(anno[, .N, by = "SUBJID"]$N, main = "Samples by subject", xlab = "# samples", xlim = c(0, 100), breaks = 100)
anno[, .N, by = "SUBJID"][N>100,]
anno[SUBJID == "K-562", SMTSD][1]
#' probably so many replicates due to big mutation rate of "common ancestor" cell

