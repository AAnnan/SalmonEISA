getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
txFile <- "rawcounts/GROseq_RNAseq2019/rawcounts_transcript.txt"
geFile <- "rawcounts/GROseq_RNAseq2019/rawcounts_gene.txt"

gro_txFile <- "rawcounts/GROseq_RNAseq2019/rawcounts_transcript_gro.txt"
gro_geFile <- "rawcounts/GROseq_RNAseq2019/rawcounts_gene_gro.txt"
ctrlRNAi <- TRUE
chrm = "all"

prepoc_cnt <- function(geFile, txFile, gene_table) {
  # read in count tables
  cntIn <- get_cnt(geFile)
  # aggregate the read counts from transcripts to genes
  cntEx <- get_cnt(txFile)
  
  genes.in.both <- intersect(rownames(cntEx),rownames(cntIn))
  cntIn <- get_names(cntIn, genes.in.both, gene_table, chrom=chrm)
  cntEx <- get_names(cntEx, genes.in.both, gene_table, chrom=chrm)
  #all(rownames(cntIn)==rownames(cntEx))
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  cntEx.norm <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx))) # normalize samples to avearge sequencing depth for exons
  cntIn.norm <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn))) # normalize samples to avearge sequencing depth for introns
  
  #genes.sel <- rowMeans(log2(cntEx.norm+8))>=4.321928 & rowMeans(log2(cntIn.norm+8))>=4.321928 #20 (12)
  #genes.sel <- rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5 #32(24)
  genes.sel <- rowSums(cntIn.norm != 0)>=(ncol(cntIn)-1) & rowSums(cntEx.norm != 0)>=(ncol(cntEx)-1) & rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5 #32(24) 
  
  # Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
  cntEx <- cntEx[genes.sel,]
  cntIn <- cntIn[genes.sel,]
  return(list(cntEx,cntIn))
}

cnt2 <- prepoc_cnt(geFile, txFile, gene_table)
cntEx <- data.frame(cnt2[1])
cntIn <- data.frame(cnt2[2])

# read in count tables
gro_cntIn <- get_cnt(gro_geFile)
# aggregate the read counts from transcripts to genes
gro_cntEx <- get_cnt(gro_txFile)

if (ctrlRNAi) {
  gro_cntEx[,c(1,2)] <- NULL
  gro_cntIn[,c(1,2)] <- NULL
} else {
  gro_cntEx[,c(3,4)] <- NULL
  gro_cntIn[,c(3,4)] <- NULL
}

gro_genes.in.both <- intersect(rownames(gro_cntEx),rownames(gro_cntIn))
gro_cntIn <- get_names(gro_cntIn, gro_genes.in.both, gene_table, chrom=chrm)
gro_cntEx <- get_names(gro_cntEx, gro_genes.in.both, gene_table, chrom=chrm)
all(rownames(gro_cntIn)==rownames(gro_cntEx))

gro_cnt <- gro_cntEx + gro_cntIn
# find genes with sufficient exonic and intronic counts (genes.sel)
gro_cnt.norm <- as.data.frame(t(mean(colSums(gro_cnt))*t(gro_cnt)/colSums(gro_cnt))) # normalize samples to avearge sequencing depth for exons

#gro_genes.sel <- rowMeans(log2(gro_cnt.norm+8))>=4.321928 #20 (12)
gro_genes.sel <- rowSums(gro_cnt.norm != 0)>=(ncol(gro_cnt.norm)-1) & rowMeans(log2(gro_cnt.norm+8))>=5#32(24) 

#gro_genes.sel <- rowMeans(log2(gro_cnt.norm+8))>=5 #32 (24)

# Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
gro_cnt <- gro_cnt[gro_genes.sel,]
gro_cnt.norm <- gro_cnt.norm[gro_genes.sel,]

# edgeR workflow
#conditions <- c("366","366","366","366","382","382","382","382")
conditions <- c("WT","WT","WT","WT","WT","sdc2RNAi","sdc2RNAi","sdc2RNAi")
factorCondition <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
group <- c(1,1,1,1,1,2,2,2)

gro_conditions <- c("WT", "WT","sdc2RNAi","sdc2RNAi")
gro_factorCondition <- factor(gro_conditions, levels=unique(gro_conditions)) # define experimental factor 'conditions'
gro_group <- c(1,1,2,2)

edgeR_flow <- function(cntEx, cntIn, groupE, factorConditionE) {
  ##Exons
  yEx <- DGEList(counts=cntEx, genes=rownames(cntEx), group=groupE) # define DGEList object
  yEx <- calcNormFactors(yEx) # determine normalization factors
  designEx <- model.matrix(~ factorConditionE) # design matrix
  rownames(designEx) <- colnames(cntEx)
  yEx <- estimateDisp(yEx, designEx) # estimate dispersion
  fitEx <- glmFit(yEx, designEx) # fit generalized linear model
  lrtEx <- glmLRT(fitEx) # calculate likelihood-ratio between full and reduced models
  ttEx <- topTags(lrtEx, n=nrow(yEx)) #final table with significance level for each gene 
  
  ##Introns
  yIn <- DGEList(counts=cntIn, genes=rownames(cntIn), group=groupE) # define DGEList object
  yIn <- calcNormFactors(yIn) # determine normalization factors
  designIn <- model.matrix(~ factorConditionE) # design matrix
  rownames(designIn) <- colnames(cntIn)
  yIn <- estimateDisp(yIn, designIn) # estimate dispersion
  fitIn <- glmFit(yIn, designIn) # fit generalized linear model
  lrtIn <- glmLRT(fitIn) # calculate likelihood-ratio between full and reduced models
  ttIn <- topTags(lrtIn, n=nrow(yIn)) #final table with significance level for each gene 

  return(list(ttEx,ttIn))
}

tt2 <- edgeR_flow(cntEx, cntIn, group, factorCondition)
ttEx <- data.frame(tt2[1])
ttIn <- data.frame(tt2[2])

##GRO
gro_y <- DGEList(counts=gro_cnt, genes=rownames(gro_cnt), group=gro_group) # define DGEList object
gro_y <- calcNormFactors(gro_y) # determine normalization factors
gro_design <- model.matrix(~ gro_factorCondition) # design matrix
rownames(gro_design) <- colnames(gro_cnt)
gro_y <- estimateDisp(gro_y, gro_design) # estimate dispersion
gro_fit <- glmFit(gro_y, gro_design) # fit generalized linear model
gro_lrt <- glmLRT(gro_fit) # calculate likelihood-ratio between full and reduced models
gro_tt <- topTags(gro_lrt, n=nrow(gro_y)) #final table with significance level for each gene 


tt.df.in <- data.frame(ttIn)
tt.df.ex <- data.frame(ttEx)
delta.cnt.edgeR <- data.frame(row.names = rownames(tt.df.ex[order(row.names(tt.df.ex)),]),dExon=tt.df.ex$logFC[order(row.names(tt.df.ex))], dIntron=tt.df.in$logFC[order(row.names(tt.df.in))]) 
forgro_delta.cnt.edgeR <- data.frame(row.names = rownames(tt.df.ex[order(row.names(tt.df.ex)),]),dIntron_WT_sdc2RNAi=tt.df.in$logFC[order(row.names(tt.df.in))]) 

gro.tt.df <- data.frame(gro_tt)
gro.delta.cnt.edgeR <- data.frame(row.names = rownames(gro.tt.df[order(row.names(gro.tt.df)),]),GROseq_WT_sdc2RNAi=gro.tt.df$logFC[order(row.names(gro.tt.df))]) 

gro_scatter(forgro_delta.cnt.edgeR,gro.delta.cnt.edgeR)


gro_delta <- data.frame(names = rownames(gro_cnt.norm),GROseq_N2_sdc2RNAi=log2(gro_cnt.norm$GROseq_sdc2/gro_cnt.norm$N2_GROseq_ctrl))
gro_delta <- gro_delta[is.finite(gro_delta$GROseq_N2_sdc2RNAi),]
rownames(gro_delta) <- gro_delta$names
gro_delta$names <- NULL
gro_scatter(forgro_delta.cnt.edgeR,gro_delta)



## Select genes with significant delta Intron/Exon (False Discovery rate inferior than 5%)
ttEx$table$PValue <- -log(p.adjust(ttEx$table$PValue, method = "BH"))
ttIn$table$PValue <- -log(p.adjust(ttIn$table$PValue, method = "BH"))
signiEx <- ttEx$table[ttEx$table$PValue>2.995732,]
signiIn <- ttIn$table[ttIn$table$PValue>2.995732,]
both_signi <- intersect(rownames(signiIn),rownames(signiEx))
delta.cnt.edgeR.signi <- delta.cnt.edgeR[rownames(delta.cnt.edgeR) %in% both_signi,]

##Plots

delta.cnt.edgeR_X <- sel_X(delta.cnt.edgeR,gene_table)
scatter_deltas_X(delta.cnt.edgeR,delta.cnt.edgeR_X)
#scatter_deltas_s3(delta.cnt.edgeR,signiEx,signiIn)

#Looking for specific genes
delta.cnt.signi[delta.cnt.signi$dIntron< -2 & delta.cnt.signi$dExon> -1 ,]
delta.cnt.edgeR[delta.cnt.edgeR$dIntron>6 | delta.cnt.edgeR$dExon>5 ,]
delta.cnt.signi[rownames(delta.cnt.signi) == "sep-1",]
cntEx.norm[rownames(cntEx.norm) == "sep-1",]
cntIn.norm[rownames(cntIn.norm) == "sep-1",]

