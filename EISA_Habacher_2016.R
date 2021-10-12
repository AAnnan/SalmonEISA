getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")
#setwd("~/Documents/Project_IZB/Papers/EISA/Habacher 2016 Data/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(edgeR)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "/Users/aa/Documents/Project_IZB/biodata/rnaseq/wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
conditions <- c("N2_mockRNAi","N2_mockRNAi","N2_rege1RNAi","N2_rege1RNAi")#,"glp1_mockRNAi","glp1_mockRNAi","glp1_rege1","glp1_rege1")
#ofInt <- "dIdE_detectedGenes.txt"
fromScratch=T
#Functions
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
scatter_deltas_H <- function(delta.cnt,delta.cnt.red) {
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  
  `%notin%` <- Negate(`%in%`)
  ets4 <- delta.cnt[(rownames(delta.cnt)=="ets-4"),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.red)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=1, size=0.7) +
    geom_point(data=delta.cnt.red, mapping=aes(x=dIntron, y=dExon, color = "Putative post-tc\nregulated genes"), alpha=1, size=0.7) +
    geom_point(data=ets4, mapping=aes(x=dIntron, y=dExon, color = "ets-4"), alpha=1, size=0.7) +
    ggtitle(paste0('R = ',corEvI)) +
    scale_colour_manual(name = NULL, values = c("Putative post-tc\nregulated genes" = "red", "ets-4" = "green3")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))
    
}
scatter_deltas_EdgeR <- function(ttIn,ttEx) {
  delta.cnt <- data.frame(dExon=tt.df.ex$logFC, dIntron=tt.df.in$logFC, row.names = rownames(tt.df.ex))
  
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  
  `%notin%` <- Negate(`%in%`)
  ets4 <- delta.cnt[(rownames(delta.cnt)=="ets-4"),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.red)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=1, size=0.7) +
    geom_point(data=delta.cnt.red, mapping=aes(x=dIntron, y=dExon, color = "Putative post-tc\nregulated genes"), alpha=1, size=0.7) +
    geom_point(data=ets4, mapping=aes(x=dIntron, y=dExon, color = "ets-4"), alpha=1, size=0.7) +
    ggtitle(paste0('R = ',corEvI)) +
    scale_colour_manual(name = NULL, values = c("Putative post-tc\nregulated genes" = "red", "ets-4" = "green3")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))
    
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
  
  cntExRawS <- cntEx
  cntInRawS <- cntIn
  
  # Normalize by library size
  cntEx <- as.data.frame(t(mean(colSums(cntEx))*t(cntEx)/colSums(cntEx))) # normalize samples to avearge sequencing depth for exons
  cntIn <- as.data.frame(t(mean(colSums(cntIn))*t(cntIn)/colSums(cntIn))) # normalize samples to avearge sequencing depth for introns
  
  # Add 8 pseudo counts
  cntEx <- log2(cntEx + 8)
  cntIn <- log2(cntIn + 8)
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  #genes.sel <- rowMeans(cntEx)>=4.321928 & rowMeans(cntIn)>=4.321928 #20 (12)
  genes.sel <- rowMeans(cntEx)>=6 & rowMeans(cntIn)>=6 #64 (56)
  
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
  
  cntExRawH <- (2^cntExD) - 8
  cntInRawH <- (2^cntInD) - 8
  
  # find genes with sufficient exonic and intronic counts (genes.sel)
  genes.sel <- rowMeans(cntExD)>=4.321928 & rowMeans(cntInD)>=4.321928 #20 (12)
  #genes.sel <- rowMeans(cntEx)>=5.247928 & rowMeans(cntIn)>=5.247928 #38(30)
  
  # fusion int/ex df
  cnt <- cbind(Ex=cntExD[genes.sel,], In=cntInD[genes.sel,])
}

# edgeR workflow
factorCondition <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
group <- c(1,1,2,2)

if (fromScratch==FALSE) {
  cntExRaw <- cntExRawH
  cntInRaw <- cntInRawH
} else if (fromScratch==TRUE){
  cntExRaw <- cntExRawS
  cntInRaw <- cntInRawS
}

##Exons
yEx <- DGEList(counts=cntExRaw, genes=rownames(cntExRaw), group=group) # define DGEList object
yEx <- calcNormFactors(yEx) # determine normalization factors
designEx <- model.matrix(~ factorCondition) # design matrix
rownames(designEx) <- colnames(cntExRaw)
yEx <- estimateDisp(yEx, designEx) # estimate dispersion
fitEx <- glmFit(yEx, designEx) # fit generalized linear model
lrtEx <- glmLRT(fitEx) # calculate likelihood-ratio between full and reduced models
ttEx <- topTags(lrtEx, n=nrow(yEx)) #final table with significance level for each gene 
head(ttEx$table)

##Introns
yIn <- DGEList(counts=cntInRaw, genes=rownames(cntInRaw), group=group) # define DGEList object
yIn <- calcNormFactors(yIn) # determine normalization factors
designIn <- model.matrix(~ factorCondition) # design matrix
rownames(designIn) <- colnames(cntInRaw)
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
delta.cnt <- get_deltas_H(cnt)

# Select putative REGE-1 targets that predominantly change (more than 2-fold) in the exonic (mature mRNA)
# but not the intronic (nascent RNA) reads
redHabach <-intersect(rownames(ttEx$table[ttEx$table$logFC>1.1,]),rownames(ttIn$table[abs(ttIn$table$logFC)<1,]))

#rg <- read.csv("/Users/aa/Documents/GitHub/SalmonEISA/red_genes.csv", header = FALSE)
#redG <- intersect(rg$V1, both_signi)

#delta.cnt.red <- delta.cnt[rownames(delta.cnt) %in% rg$V1,] #Select red genes
delta.cnt.red <- delta.cnt[rownames(delta.cnt) %in% redHabach,] #Select red genes
##Plots
scatter_deltas_H(delta.cnt,delta.cnt.red)

delta.cnt.black <- delta.cnt[rownames(delta.cnt) %in% both_signi,] #Select red genes
scatter_deltas_H2 <- function(delta.cnt,delta.cnt.red,delta.cnt.black) {
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  
  `%notin%` <- Negate(`%in%`)
  ets4 <- delta.cnt[(rownames(delta.cnt)=="ets-4"),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.red)),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.black)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=0.5, size=1) +
    geom_point(data=delta.cnt.black, mapping=aes(x=dIntron, y=dExon, color = "FDR<0.05\ndIntron or dExon\n"), alpha=0.5, size=1) +
    geom_point(data=delta.cnt.red, mapping=aes(x=dIntron, y=dExon, color = "Putative post-tc\nregulated genes\n"), alpha=0.5, size=1) +
    geom_point(data=ets4, mapping=aes(x=dIntron, y=dExon, color = "ets-4"), alpha=0.5, size=1.2) +
    ggtitle(paste0('R = ',corEvI)) +
    scale_colour_manual(name = NULL, values = c("FDR<0.05\ndIntron or dExon\n"="blue","Putative post-tc\nregulated genes\n" = "red", "ets-4" = "green3")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))
    
}
scatter_deltas_s3 <- function(delta.cnt,signiEx,signiIn) {
  `%notin%` <- Negate(`%in%`)
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  both_signi <- intersect(rownames(signiIn),rownames(signiEx))
  signiIn <- signiIn[rownames(signiIn) %notin% both_signi,]
  signiEx <- signiEx[rownames(signiEx) %notin% both_signi,]
  
  delta.cnt.Sin <- delta.cnt[rownames(delta.cnt) %in% rownames(signiIn),]
  delta.cnt.Sex <- delta.cnt[rownames(delta.cnt) %in% rownames(signiEx),]
  delta.cnt.signi <- delta.cnt[rownames(delta.cnt) %in% both_signi,]
  
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.Sin)),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.Sex)),]
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.signi)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=0.5, size=1) +
    geom_point(data=delta.cnt.Sex, mapping=aes(x=dIntron, y=dExon, color = "FDR<0.05\ndExon only\n"), alpha=0.5, size=1) +
    geom_point(data=delta.cnt.Sin, mapping=aes(x=dIntron, y=dExon, color = "FDR<0.05\ndIntron only\n"), alpha=0.5, size=1) +
    geom_point(data=delta.cnt.signi, mapping=aes(x=dIntron, y=dExon, color = "FDR<0.05\ndIntron and dExon\n"), alpha=0.7, size=1) +
    ggtitle(paste0('R = ',corEvI)) +
    scale_colour_manual(name = NULL, values = c("FDR<0.05\ndExon only\n" = "yellow3", "FDR<0.05\ndIntron only\n" = "green3","FDR<0.05\ndIntron and dExon\n"="red")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))
    
}
#scatter_deltas_s3(delta.cnt,signiEx,signiIn)
##Plots
scatter_deltas_H2(delta.cnt,delta.cnt.red, delta.cnt.black)

tt.df.in <- data.frame(ttIn)
tt.df.ex <- data.frame(ttEx)

tt.df.in <- tt.df.in[ order(row.names(tt.df.in)), ]
tt.df.ex <- tt.df.ex[ order(row.names(tt.df.ex)), ]

if (fromScratch==FALSE) {
  genes.sel2 <- tt.df.ex$logCPM>=4 & tt.df.in$logCPM>=4 #20 (12)
} else if (fromScratch==TRUE){
  genes.sel2 <- tt.df.ex$logCPM>=6 & tt.df.in$logCPM>=6 #64 (56)
}

tt.df.in <- tt.df.in[genes.sel2,]
tt.df.ex <- tt.df.ex[genes.sel2,]

scatter_deltas_EdgeR(tt.df.in,tt.df.ex)
delta.cnt.EdgeR <- data.frame(dExon=tt.df.ex$logFC, dIntron=tt.df.in$logFC, row.names = rownames(tt.df.ex))
scatter_deltas_s3(delta.cnt.EdgeR,signiEx,signiIn)

delta.cnt.edgeR.signi[delta.cnt.edgeR.signi$dIntron>2 & delta.cnt.edgeR.signi$dExon>2 ,]
