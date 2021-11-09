getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

##### Files (counts, geneTable) and parameters (EdgeR)
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
txFile <- "rawcounts/REST_ALL/WT2020_vs_scc1_2020_rawcounts_transcript.txt"
geFile <- "rawcounts/REST_ALL/WT2020_vs_scc1_2020_rawcounts_gene.txt"

conditions <- c("WT","WT","WT","WT","dpy26","dpy26","dpy26","dpy26")
#######



# read in the "genes" counts as intronic
cntIn <- read.delim(geFile, row.names = 1)
# aggregate the transcripts counts by genes, as exonic
cntEx <- get_cnt(txFile)

#list of WORMBASE IDs present in both introns and exons raw counts
genes.in.both <- intersect(rownames(cntEx),rownames(cntIn))
# Get rid of genes existing only in 1 count table
cntIn <- cntIn[rownames(cntIn) %in% genes.in.both,]
cntEx <- cntEx[rownames(cntEx) %in% genes.in.both,]

# Change Gene_IDs to gene names
cntIn <- get_names(cntIn, gene_table)
cntEx <- get_names(cntEx, gene_table)

# Normalize exonic and intronic counts to av. sequencing depth and filter by abundance
# find genes with sufficient exonic and intronic counts (genes.sel)
cntEx.norm <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx)))
cntIn.norm <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn)))

genes.sel <- rowSums(cntIn.norm != 0)>=7 & rowSums(cntEx.norm != 0)>=7 & rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5 #32(24) 
#4.321928 #20 (12)

# Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
cntEx <- cntEx[genes.sel,]
cntIn <- cntIn[genes.sel,]

# edgeR workflow
ttEx <- edgeR_flow(cntEx,conditions)
ttIn <- edgeR_flow(cntIn,conditions)

delta.cnt <- data.frame(row.names = rownames(ttEx),
                        dExon=ttEx$logFC,
                        dIntron=ttIn$logFC)


delta.cnt.X <- sel_chrom(delta.cnt, gene_table, "X")
scatter_deltas_X(delta.cnt, delta.cnt.X, conditions)

#Looking for outlier genes on the scatter plot
print(delta.cnt[delta.cnt$dIntron <= -1 & delta.cnt$dExon >= -5 ,])
