getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
txFile <- "rawcounts/bolaji1612_rnaseq/ALL.txt"
conditions <- c("TIR1m","TIR1m","TIR1m","TIR1p","TIR1p","TIR1p",
                "TOP1m","TOP1m","TOP1m","TOP1p","TOP1p","TOP1p",
                "TOP2m","TOP2m","TOP2m","TOP2p","TOP2p","TOP2p",
                "N2","N2","AMA","AMA","FLAVO","FLAVO")

# aggregate the transcripts counts by genes, as exonic
cntEx <- get_cnt(txFile)

# Change Gene_IDs to gene names
cntEx <- get_names(cntEx, gene_table)

# Normalize exonic and intronic counts to av. sequencing depth and filter by abundance
# find genes with sufficient exonic and intronic counts (genes.sel)
cntEx.norm <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx)))

genes.sel <- rowSums(cntEx.norm != 0)>=(ncol(cntEx)) & rowMeans(log2(cntEx.norm+8))>=4.584963 #24 (16) #32(24) 
#4.321928 #20 (12)

# Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
cnt <- cntEx[genes.sel,]

COND_EXP <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
#group <- c(1,1,2,2,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9)
#y <- DGEList(counts=cnt, genes=rownames(cnt), group=group) # define DGEList object
y <- DGEList(counts=cnt, genes=rownames(cnt))
y <- calcNormFactors(y) # determine normalization factors

design <- model.matrix(~ COND_EXP) # design matrix
rownames(design) <- colnames(cnt)
y <- estimateDisp(y, design) # estimate dispersion
y$design
plotMDS(y)