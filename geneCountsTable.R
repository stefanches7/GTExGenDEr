setwd("~/projects/GTExAna/")
library(data.table)
dt <- fread(gzfile("zcat < data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"))
