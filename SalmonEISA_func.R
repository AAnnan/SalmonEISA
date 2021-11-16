### Functions used in SalmonEISA.R

#### Imports ####

#' Extracts a read count per gene df from Salmon's output via 
#' the files output by file_prep_Salmon.sh
#' 
#' @param path_rc path to a tsv listing read counts per exons (line number) per genes (line)
#' @return A large matrix with read counts summed by genes (rownames).
#' @examples
#' cntEx <- get_cnt("rawcounts_intronic.txt")
get_cnt <- function(path_rc) {
  rc <- read.delim(path_rc)
  colnames(rc)[1] <- "Gene_IDs"
  cnts <- aggregate(.~Gene_IDs, data=rc, FUN=sum)
  
  rownames(cnts) <- cnts[,1]
  cnts[,1] <- NULL
  
  return(cnts)
}

#' Replace Gene IDs by gene names
#' 
#' @param cnt A df with WORMBASE IDs as rownames
#' @param path_table path to a look up table with gene names, WORMBASE IDs and chromosomes
#' @param chrom select a chromosome (ex:X) or all
#' @return A df with Gene names as rownames
#' @examples
#' cntEx <- get_names(cntEx, gene_table, chrom="all")
get_names <- function(cnt, path_table) {
  
  # Import Gene name/chromosome look up table
  gene.table <- read.delim2(path_table, row.names = 1)
  #Build gene name list from gene table, exclude genes not found (NAs)
  ge_list <- gene.table$Gene.name[match(rownames(cnt), rownames(gene.table))] 
  #only include "not NA rows" in cnt
  cnt <- cnt[!is.na(ge_list),]
  #replace rownames Gene_IDs by names
  rownames(cnt) <- na.omit(ge_list)
  
  return(cnt)
}


#### Manipulations ####
# edgeR workflow
edgeR_flow <- function(cnt, conditions) {
  COND_EXP <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
  ##Exons
  y <- DGEList(counts=cnt, genes=rownames(cnt)) # define DGEList object
  y <- calcNormFactors(y) # determine normalization factors
  design <- model.matrix(~ COND_EXP) # design matrix
  rownames(design) <- colnames(cnt)
  y <- estimateDisp(y, design) # estimate dispersion
  fit <- glmFit(y, design) # fit generalized linear model
  lrt <- glmLRT(fit) # calculate likelihood-ratio between full and reduced models
  tt <- data.frame(row.names = 1, topTags(lrt, n=nrow(y), adjust.method = "BH", sort.by = "none")) #final table with significance level for each gene 
  
  return(tt[order(row.names(tt)),])
}

sel_chrom <- function(df, path_table, chrom="X") {
  gene.table <- read.delim2(row.names = 1, path_table)
  gene.table <- gene.table[gene.table$Chromosome.scaffold.name==chrom,]
  
  ge <- gene.table$Gene.name[match(rownames(df), gene.table$Gene.name)]
  df <- df[!is.na(ge),]
  rownames(df) <- na.omit(ge)
  
  return(df)
}

#### Plots ####

#' Build scatter plot of DeltaIntron/Exon 
#' 
#' @param delta.cnt An eisa df with gene names as rownames and at least 2 columns:
#' $dIntron and $dExon
#' @param delta.cnt_X A subset eisa df with gene names as rownames and at least 2 columns:
#' $dIntron and $dExon
#' @return A scatter plot of DeltaIntron/Exon with Pearson correlation coeff
#' @examples
#' scatter_deltas(delta.cnt,delta.cnt.signi)
scatter_deltas_X <- function(delta.cnt,delta.cnt_X,conditions) {
  `%notin%` <- Negate(`%in%`)
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  corEvI_X <- round(cor(delta.cnt_X$dIntron,delta.cnt_X$dExon), 3)
  
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt_X)),]
  exp_cond <- paste(rev(unique(conditions)), collapse=' - ' )
  
  plt <- ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon, color = "Autosomal\ngenes"), alpha=0.5, size=1) +
    geom_point(data=delta.cnt_X, mapping=aes(x=dIntron, y=dExon, color = "X genes"), alpha=0.5, size=1) +
    ggtitle(paste0(exp_cond,'\n\nR = ',corEvI), subtitle = paste0('R = ',corEvI_X)) +
    scale_colour_manual(name = NULL, values = c("X genes"="red","Autosomal\ngenes"="black")) +
    labs(x= "\u0394Intron",    # works fine
         y= "\u0394Exon") +  # works fine
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.minor = element_blank()) +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=45, b = -40)),
          plot.subtitle=element_text(size=12, face="italic", color="red", margin = margin(t=45, b = -40)))
  
  t <- round(max(layer_scales(plt)$y$range$range[2],layer_scales(plt)$x$range$range[2]),1)
  s <- round(min(layer_scales(plt)$y$range$range[1],layer_scales(plt)$x$range$range[1]),1)
  
  plt <- plt + coord_equal(xlim=c(s,t), ylim=c(s,t)) + 
    geom_segment(aes(x = s, y = s,xend = t, yend = t),size = 0.2) +
    theme(plot.margin = unit(c(-1,-2,0.3,-2), "cm"))
  
  plt
  #ggsave(paste0(exp_cond,".svg"), plot = plt)
}

### ALL ###

filtering <- function(cnt){
  # Normalize exonic and intronic counts to av. sequencing depth and filter by abundance
  # find genes with sufficient exonic and intronic counts (genes.sel)
  cnt.norm <- as.data.frame(t(mean(colSums(cnt))*t(cnt)/colSums(cnt)))
  
  genes.sel <- rowSums(cnt.norm != 0)>=(ncol(cnt)-1) & rowMeans(log2(cnt.norm+8))>=5 #32(24) 
  #4.321928 #20 (12)
  
  # Keep counts only of genes with sufficient exonic and intronic counts (genes.sel)
  return(cnt[genes.sel,])
}
edgeR_flow <- function(cnt, conditions) {
  COND_EXP <- factor(conditions, levels=unique(conditions)) # define experimental factor 'conditions'
  ##Exons
  y <- DGEList(counts=cnt, genes=rownames(cnt)) # define DGEList object
  y <- calcNormFactors(y) # determine normalization factors
  design <- model.matrix(~ COND_EXP) # design matrix
  rownames(design) <- colnames(cnt)
  y <- estimateDisp(y, design) # estimate dispersion
  fit <- glmFit(y, design) # fit generalized linear model
  lrt <- glmLRT(fit) # calculate likelihood-ratio between full and reduced models
  tt <- data.frame(row.names = 1, topTags(lrt, n=nrow(y), adjust.method = "BH", sort.by = "none")) #final table with significance level for each gene 
  
  tt <- tt[order(row.names(tt)),]
  
  dfTT <- data.frame(logFC=tt$logFC, log_Adj_Pval=-log(tt$PValue))
  rownames(dfTT)<-rownames(tt)
  return(dfTT)
}
sel_rnaseqX <- function(rnaseq,gene.table) {
  gene.table <- read.delim2(gene.table)
  rownames(gene.table) <- gene.table[,1]
  gene.table <- gene.table[gene.table$Chromosome.scaffold.name=='X',]
  ge <- gene.table$Gene.name[match(rownames(rnaseq), gene.table$Gene.name)]
  rnaseq <- rnaseq[!is.na(ge),]
  rownames(rnaseq) <- na.omit(ge)
  return(rnaseq)
}
volc_rnaseq <- function(rnaseq) {
  rnaseq<-rnaseq[rnaseq$log_Adj_Pval<200,]
  rnaseqsig <- rnaseq[-rnaseq$log_Adj_Pval<log(0.05),]
  ggplot() + 
    geom_point(data=rnaseq, mapping=aes(x=logFC, y=log_Adj_Pval),size=0.5, alpha=0.5) +
    geom_point(data=rnaseqsig, mapping=aes(x=logFC, y=log_Adj_Pval, color = "adjPval<0.05"),size=0.5, alpha=0.75) +
    ggtitle('sdc2RNAi - WT') +
    scale_colour_manual(name = NULL, values = c("adjPval<0.05" = "red")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))# +lims(x = c(-2, 2), y = c(-2, 2))
}
scatter_RNAs <- function(rnaseqA,rnaseqB) {
  rnaseqA <- rnaseqA[rownames(rnaseqA) %in% rownames(rnaseqB),]
  rnaseqB <- rnaseqB[rownames(rnaseqB) %in% rownames(rnaseqA),]
  rnaseqA <- rnaseqA[order(row.names(rnaseqA)), ]
  rnaseqB <- rnaseqB[order(row.names(rnaseqB)), ]
  
  scatR <- data.frame(logFCA=rnaseqA$logFC, logFCB=rnaseqB$logFC)
  
  corEvI <- round(cor(scatR$logFCA,scatR$logFCB), 3)
  
  ggplot() + 
    geom_point(data=scatR, mapping=aes(x=logFCA, y=logFCB), alpha=0.5) +
    ggtitle(paste0('R = ',corEvI)) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)))# +lims(x = c(-2, 2), y = c(-2, 2))
}
