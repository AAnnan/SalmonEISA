getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
txFile <- "rawcounts/24_09_GT/rawcounts_transcript.txt"
geFile <- "rawcounts/24_09_GT/rawcounts_gene.txt"
#conditions <- c("366","366","366","366","382","382","382","382")
conditions <- c("WT","WT","WT","WT","dpy26cs","dpy26cs","dpy26cs","dpy26cs")

# read in count tables
#cntGe <- get_cnt(geFile)
cntIn <- read.delim(geFile, row.names = 1)
# aggregate the read counts from transcripts to genes
cntEx <- get_cnt(txFile)

# Get rid of genes existing only in 1 count table
# Change Gene_IDs to gene names
#genes.in.both <- intersect(rownames(cntEx),rownames(cntGe))
#cntGe <- cntGe[rownames(cntGe) %in% genes.in.both,]
#cntEx <- cntEx[rownames(cntEx) %in% genes.in.both,]

#cntGe <- cntGe[ order(row.names(cntGe)), ]
#cntEx <- cntEx[ order(row.names(cntEx)), ]

#cntIn <- cntGe - cntEx
#cntIn[cntIn < 0] <- NA
#cntIn <- na.omit(cntIn)

genes.in.both <- intersect(rownames(cntEx),rownames(cntIn))
cntIn <- get_names(cntIn, genes.in.both, gene_table, chrom="all")
cntEx <- get_names(cntEx, genes.in.both, gene_table, chrom="all")

# find genes with sufficient exonic and intronic counts (genes.sel)
cntEx.norm <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx))) # normalize samples to avearge sequencing depth for exons
cntIn.norm <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn))) # normalize samples to avearge sequencing depth for introns

#genes.sel <- rowMeans(log2(cntEx.norm+8))>=4.321928 & rowMeans(log2(cntIn.norm+8))>=4.321928 #20 (12)
genes.sel <- rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5 #32(24)

# combine exon and intron raw counts for edgeR containing genes with sufficient counts in both exonic and intronic levels
cnt.norm <- cbind(Ex=cntEx.norm[genes.sel,], In=cntIn.norm[genes.sel,])
#cnt.norm.norm <- as.data.frame(t(mean(colSums(cnt.norm))*t(cnt.norm)/colSums(cnt.norm)))
#cnt <- cbind(Ex=cntEx[genes.sel,], In=cntIn[genes.sel,])

# Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
cntEx <- cntEx[genes.sel,]
cntIn <- cntIn[genes.sel,]

# edgeR workflow
factorCondition <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
group <- c(1,1,1,1,2,2,2,2)

##Exons
yEx <- DGEList(counts=cntEx, genes=rownames(cntEx), group=group) # define DGEList object
yEx <- calcNormFactors(yEx) # determine normalization factors
designEx <- model.matrix(~ factorCondition) # design matrix
rownames(designEx) <- colnames(cntEx)
yEx <- estimateDisp(yEx, designEx) # estimate dispersion
fitEx <- glmFit(yEx, designEx) # fit generalized linear model
lrtEx <- glmLRT(fitEx) # calculate likelihood-ratio between full and reduced models
ttEx <- topTags(lrtEx, n=nrow(yEx)) #final table with significance level for each gene 
head(ttEx$table)

##Introns
yIn <- DGEList(counts=cntIn, genes=rownames(cntIn), group=group) # define DGEList object
yIn <- calcNormFactors(yIn) # determine normalization factors
designIn <- model.matrix(~ factorCondition) # design matrix
rownames(designIn) <- colnames(cntIn)
yIn <- estimateDisp(yIn, designIn) # estimate dispersion
fitIn <- glmFit(yIn, designIn) # fit generalized linear model
lrtIn <- glmLRT(fitIn) # calculate likelihood-ratio between full and reduced models
ttIn <- topTags(lrtIn, n=nrow(yIn)) #final table with significance level for each gene 
head(ttIn$table)

## Select genes with significant delta Intron/Exon (False Discovery rate inferior than 5%)
signiEx <- ttEx$table[ttEx$table$FDR<0.05,]
signiIn <- ttIn$table[ttIn$table$FDR<0.05,]
both_signi <- union(rownames(signiIn),rownames(signiEx))

## Average over replicates, build delta Intron/Exon with error bars (mean+-sd)
delta.cnt <- get_deltas(cnt.norm,logBefore=TRUE)
delta.cnt.signi <- delta.cnt[rownames(delta.cnt) %in% both_signi,] #Select significant genes

##Plots
#plot_col_deltas(delta.cnt)
#plot_col_deltas(delta.cnt.signi)

scatter_deltas(delta.cnt,delta.cnt.signi)

## Average over replicates, build mean Intron/Exon counts with error bars (mean+-sd)
#mean.cnt <- get_means(cnt.norm)
#mean.cnt.signi <- mean.cnt[rownames(mean.cnt) %in% both_signi,] #Select significant genes
#plot_col_means(mean.cnt)
#plot_col_means(mean.cnt.signi)
