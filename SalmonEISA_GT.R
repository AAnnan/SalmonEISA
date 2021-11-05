getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/rawcounts/ALL_GT/")

library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "../../wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"

Salmoneisa <- function(gene_table,geFile,txFile,conditions,group,cond) {
  # read in count tables
  cntIn <- get_cnt(geFile)
  cntEx <- get_cnt(txFile)
  
  genes.in.both <- intersect(rownames(cntEx),rownames(cntIn))
  cntIn <- get_names(cntIn, genes.in.both, gene_table, chrom="all")
  cntEx <- get_names(cntEx, genes.in.both, gene_table, chrom="all")
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  cntEx.norm <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx))) # normalize samples to avearge sequencing depth for exons
  cntIn.norm <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn))) # normalize samples to avearge sequencing depth for introns
  
  #genes.sel <- rowMeans(log2(cntEx.norm+8))>=4.321928 & rowMeans(log2(cntIn.norm+8))>=4.321928 #20 (12)
  genes.sel <- rowMeans(log2(cntEx.norm+8))>=5 & rowMeans(log2(cntIn.norm+8))>=5 #32(24)
  
  # Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
  cntEx <- cntEx[genes.sel,]
  cntIn <- cntIn[genes.sel,]
  
  #EdgeR
  factorCondition <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
  ##Exons
  yEx <- DGEList(counts=cntEx, genes=rownames(cntEx), group=group) # define DGEList object
  yEx <- calcNormFactors(yEx) # determine normalization factors
  designEx <- model.matrix(~ factorCondition) # design matrix
  rownames(designEx) <- colnames(cntEx)
  yEx <- estimateDisp(yEx, designEx) # estimate dispersion
  fitEx <- glmFit(yEx, designEx) # fit generalized linear model
  lrtEx <- glmLRT(fitEx) # calculate likelihood-ratio between full and reduced models
  ttEx <- topTags(lrtEx, n=nrow(yEx)) #final table with significance level for each gene 
  
  ##Introns
  yIn <- DGEList(counts=cntIn, genes=rownames(cntIn), group=group) # define DGEList object
  yIn <- calcNormFactors(yIn) # determine normalization factors
  designIn <- model.matrix(~ factorCondition) # design matrix
  rownames(designIn) <- colnames(cntIn)
  yIn <- estimateDisp(yIn, designIn) # estimate dispersion
  fitIn <- glmFit(yIn, designIn) # fit generalized linear model
  lrtIn <- glmLRT(fitIn) # calculate likelihood-ratio between full and reduced models
  ttIn <- topTags(lrtIn, n=nrow(yIn)) #final table with significance level for each gene 
  
  tt.df.in <- data.frame(ttIn)
  tt.df.ex <- data.frame(ttEx)
  delta.cnt.edgeR <- data.frame(row.names = rownames(tt.df.ex[order(row.names(tt.df.ex)),]),dExon=tt.df.ex$logFC[order(row.names(tt.df.ex))], dIntron=tt.df.in$logFC[order(row.names(tt.df.in))]) 
  
  delta.cnt.edgeR_X <- sel_X(delta.cnt.edgeR,gene_table)

  #Looking for specific genes
  print(delta.cnt.edgeR[rownames(delta.cnt.edgeR)=="srw-85",])
  
  scatter_deltas_X(delta.cnt.edgeR,delta.cnt.edgeR_X,cond)
}

Salmoneisa(gene_table,"WT2020_vs_dpy26cs2020_rawcounts_gene.txt","WT2020_vs_dpy26cs2020_rawcounts_transcript.txt",c("WT","WT","WT","WT","exp","exp","exp","exp"),c(1,1,1,1,2,2,2,2),"WT2020-dpy26cs2020")
Salmoneisa(gene_table,"TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_gene.txt","TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp"),c(1,1,1,2,2,2),"TIR1sd3deg-dpy26TIR1sd3degA")
Salmoneisa(gene_table,"TIR1sd3deg_vs_TIR1sd3degA_rawcounts_gene.txt","TIR1sd3deg_vs_TIR1sd3degA_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp"),c(1,1,1,2,2,2),"TIR1sd3deg-TIR1sd3degA")
Salmoneisa(gene_table,"TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_gene.txt","TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp"),c(1,1,1,2,2,2),"TIR1wtA-dpy26TIR1sd3degA")
Salmoneisa(gene_table,"TIR1wtA_vs_TIR1sd3degA_rawcounts_gene.txt","TIR1wtA_vs_TIR1sd3degA_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp"),c(1,1,1,2,2,2),"TIR1wtA-TIR1sd3degA")


#Salmoneisa(gene_table,"WT2021aWT2020_vs_dpy26cs2020_rawcounts_gene.txt","WT2021aWT2020_vs_dpy26cs2020_rawcounts_transcript.txt",c("WT","WT","WT","WT","WT","WT","WT","exp","exp","exp","exp"),c(1,1,1,1,1,1,1,2,2,2,2),"WT2021aWT2020-dpy26cs2020")
#Salmoneisa(gene_table,"WT2021_vs_dpy26cs2020_rawcounts_gene.txt","WT2021_vs_dpy26cs2020_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp","exp"),c(1,1,1,2,2,2,2),"WT2021-dpy26cs2020")













geFiles <- list("WT2020_vs_dpy26cs2020_rawcounts_gene.txt",
                "TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",
                "TIR1sd3deg_vs_TIR1sd3degA_rawcounts_gene.txt",
                "TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",
                "TIR1wtA_vs_TIR1sd3degA_rawcounts_gene.txt",
                "WT2021aWT2020_vs_dpy26cs2020_rawcounts_gene.txt",
                "WT2021_vs_dpy26cs2020_rawcounts_gene.txt")

txFiles <- list("WT2020_vs_dpy26cs2020_rawcounts_transcript.txt",
                "TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1sd3deg_vs_TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1wtA_vs_TIR1sd3degA_rawcounts_transcript.txt",
                "WT2021aWT2020_vs_dpy26cs2020_rawcounts_transcript.txt",
                "WT2021_vs_dpy26cs2020_rawcounts_transcript.txt")

conditions <- list(c("WT","WT","WT","WT","exp","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","WT","WT","WT","WT","exp","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp","exp"))

group_s <- list(c(1,1,1,1,2,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,1,1,1,1,2,2,2,2),
                c(1,1,1,2,2,2,2))

conds <- list("WT2020-dpy26cs2020",
              "TIR1sd3deg-dpy26TIR1sd3degA",
              "TIR1sd3deg-TIR1sd3degA",
              "TIR1wtA-dpy26TIR1sd3degA",
              "TIR1wtA-TIR1sd3degA",
              "WT2021aWT2020-dpy26cs2020",
              "WT2021-dpy26cs2020")
