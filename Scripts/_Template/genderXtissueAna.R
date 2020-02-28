counts <- dt_tpm[,3:ncol(dt_tpm)]
commonTissues <- intersect(anno[SEX=="F",SMTSD], anno[SEX=="M",SMTSD])
commonSamples <- which(colnames(counts) %in% anno[SMTSD %in% commonTissues, SAMPID])
#exclude gender-specific tissues
counts <- counts[, ..commonSamples]

genders <- sapply(colnames(counts), function(x) anno[SAMPID == x, SEX])
tissues <- sapply(colnames(counts), function(x) anno[SAMPID == x, SMTSD])
tissues <- factor(tissues)
all(tissues %in% commonTissues)
design <- model.matrix(~0 + genders*tissues)

geneAnno <- fread(snakemake@input[["geneAnno"]])
#KILLZ UR PC
#y <- DGEList(counts = (counts %>% as.matrix),genes = geneAnno)
#y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
#y <- estimateGLMTrendedDisp(y, design)
#y <- estimateGLMTagwiseDisp(y, design)
#fit <- glmFit(y, design)
#lrt <- glmLRT(fit)
topTags(lrt)