source('.wBuild/wBuildParser.R')
parseWBHeader("Scripts/geneCountsTable.R")


setwd("~/projects/GTExAna/")
library(data.table)
library(ggplot2)
library(MASS)
library(magrittr)
library(edgeR)
library(tweeDEseqCountData)
library(DT)

dt_tpm <- fread(cmd = paste0("zcat < ", snakemake@input[["countsTpm"]]))

samples <- fread(snakemake@input[["sampleAnno"]])[SAMPID %in% colnames(dt_tpm)]
samples[, SUBJID:=sub('(^[^-]+-[^-]+)-(.*)$', '\\1', samples$SAMPID)] # key before the second dash
phenotypes <- fread(snakemake@input[["phenotypeAnno"]])
anno <- merge(samples, phenotypes, by = "SUBJID", all.x = T)
genderlev <- factor(anno$SEX)
levels(genderlev) <- c("M", "F")
anno[, SEX:=genderlev] 

hist(anno[, .N, by = "SUBJID"]$N, main = "Samples by subject", xlab = "# samples", xlim = c(0, 100), breaks = 100)

dt_tpm[2,5]

setdiff(tissueBySex[SEX=="F",SMTSD], tissueBySex[SEX=="M",SMTSD])

setdiff(tissueBySex[SEX=="M",SMTSD], tissueBySex[SEX=="F",SMTSD])

nrow(anno[is.na(SMTSD),]) == 0

commonTissuesAnno <- anno[SMTSD %in% intersect(tissueBySex[SEX=="F",SMTSD], tissueBySex[SEX=="M",SMTSD]), ]


tissueBySex <- anno[, .N, by = list(SMTSD, SEX)]
ggplot(tissueBySex, aes(x = SMTSD, y = N, fill = SEX)) + geom_bar(stat = "identity") +
  labs(title = "Different tissues in genders",x = "", y = "Count") + 
  theme_light() +
  theme(
    plot.margin = margin(0, 1, 1.2, 1, "cm"), 
    axis.text.x = element_text(face="plain", 
                               size=7, angle=45, hjust = 1, vjust = 0.5,
                               margin = margin(-1.7, 0, 0, 0, "cm")))