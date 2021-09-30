getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")
#setwd("~/Documents/Project_IZB/Papers/EISA/Habacher 2016 Data/")

library(ggplot2)
library(tidyr)
library(dplyr)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "/Users/aa/Documents/Project_IZB/biodata/rnaseq/wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
conditions <- c("N2_mockRNAi","N2_mockRNAi","N2_rege1RNAi","N2_rege1RNAi")#,"glp1_mockRNAi","glp1_mockRNAi","glp1_rege1","glp1_rege1")
#ofInt <- "dIdE_detectedGenes.txt"
fromScratch=FALSE

#Functions
scatter_deltas_H <- function(delta.cnt,delta.cnt.red) {
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  
  `%notin%` <- Negate(`%in%`)
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.red)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=1, size=0.7) +
    geom_point(data=delta.cnt.red, mapping=aes(x=dIntron, y=dExon, color = "ets-4"), alpha=1, size=0.7) +
    ggtitle(paste0('R = ',corEvI)) +
    scale_colour_manual(name = NULL, values = c("ets-4" = "red")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38))) +
    lims(x = c(-4.2, 8.2), y = c(-4.2, 8.2))
}
get_deltas_H <- function(cnt) {
  
  cntEx.mock <- rowMeans(cnt[,1:2])
  cntEx.rege <- rowMeans(cnt[,3:4])
  cntIn.mock <- rowMeans(cnt[,5:6])
  cntIn.rege <- rowMeans(cnt[,7:8])
  
  varEx.mock <- RowVar(cnt[,1:2])
  varEx.rege <- RowVar(cnt[,3:4])
  varIn.mock <- RowVar(cnt[,5:6])
  varIn.rege <- RowVar(cnt[,7:8])
  
  
  sd.dIntron <- sqrt(varIn.mock+varIn.rege)
  sd.dExon <- sqrt(varEx.mock+varEx.rege)
  
  dExon <- cntEx.rege-cntEx.mock
  dIntron <- cntIn.rege-cntIn.mock
  
  
  delta.cnt <- data.frame(dExon=dExon,
                          dIntron=dIntron,
                          lowI=dIntron-sd.dIntron/2,
                          highI=dIntron+sd.dIntron/2,
                          lowE=dExon-sd.dExon/2,
                          highE=dExon+sd.dExon/2)
  
  return(na.omit(delta.cnt[is.finite(rowSums(delta.cnt)),]))
  
}


if (fromScratch==TRUE) {
  
  insFile <- "rawcounts/habacher2016/rawcounts_gene.txt"
  exsFile <- "rawcounts/habacher2016/rawcounts_transcript.txt"
  
  # read in count tables and 
  # aggregate the read counts from exons/introns to genes
  cntEx <- get_cnt(exsFile)
  cntIn <- get_cnt(insFile)
  #aa <- get_cnt(insFile)
  
  #Select N2
  cntEx <- cntEx[,1:4]
  cntIn <- cntIn[,1:4]
  
  # Get rid of genes existing only in 1 count table
  # Change Gene_IDs to gene names
  genes.in.both <- intersect(rownames(cntEx),rownames(cntIn))
  
  cntEx <- get_names(cntEx, genes.in.both, gene_table, chrom="all")
  cntIn <- get_names(cntIn, genes.in.both, gene_table, chrom="all")
  
  # Normalize by library size
  cntEx <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx))) # normalize samples to avearge sequencing depth for exons
  cntIn <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn))) # normalize samples to avearge sequencing depth for introns
  
  # Add 8 pseudo counts
  cntEx <- log2(cntEx + 8)
  cntIn <- log2(cntIn + 8)
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  #genes.sel <- rowMeans(cntEx)>=4.321928 & rowMeans(cntIn)>=4.321928 #20 (12)
  genes.sel <- rowMeans(cntEx)>=6 & rowMeans(cntIn)>=6 #38(30)
  
  # fusion int/ex df
  cnt <- cbind(Ex=cntEx[genes.sel,], In=cntIn[genes.sel,])
} else if (fromScratch==FALSE) {
  #Treat read count directly from the supplementary data
  insFile <- "/Users/aa/Documents/GitHub/SalmonEISA/intronic_ets4.txt"
  exsFile <- "/Users/aa/Documents/GitHub/SalmonEISA/exonic_ets4.txt"
  
  # read in count tables and 
  # aggregate the read counts from exons/introns to genes
  cntExD <- read.delim(exsFile, row.names = 1)
  cntExD <- cntExD[ order(row.names(cntExD)), ]
  cntInD <- read.delim(insFile, row.names = 1)
  cntInD <- cntInD[ order(row.names(cntInD)), ]
  #aa <- get_cnt(insFile)
  
  #Select N2
  cntExD <- cntExD[,1:4]
  cntInD <- cntInD[,1:4]
  
  # Get rid of genes existing only in 1 count table
  # Change Gene_IDs to gene names
  genes.in.both <- intersect(rownames(cntExD),rownames(cntInD))
  
  cntExD <- get_names(cntExD, genes.in.both, gene_table, chrom="all")
  cntInD <- get_names(cntInD, genes.in.both, gene_table, chrom="all")
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  genes.sel <- rowMeans(cntExD)>=4.321928 & rowMeans(cntInD)>=4.321928 #20 (12)
  #genes.sel <- rowMeans(cntEx)>=5.247928 & rowMeans(cntIn)>=5.247928 #38(30)
  
  # fusion int/ex df
  cnt <- cbind(Ex=cntExD[genes.sel,], In=cntInD[genes.sel,])
}

## Average over replicates, build delta Intron/Exon with error bars (mean+-sd)
delta.cnt <- get_deltas_H(cnt)

# Get red genes
rg <- read.csv("/Users/aa/Documents/GitHub/SalmonEISA/red_genes.csv", header = FALSE)
delta.cnt.red <- delta.cnt[rownames(delta.cnt) %in% rg$V1,] #Select red genes

##Plots
scatter_deltas_H(delta.cnt,delta.cnt.red)



