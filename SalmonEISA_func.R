### Functions used in SalmonEISA.R

#### Stats ####
##Row operations

RowSD <- function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/dim(x)[2] - 1)
}
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
RowCov <- function(x, y){
  rowSums((x - rowMeans(x))*(y - rowMeans(y)))/(dim(x)[2] - 1)
}

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

#' Extracts a read count per gene df from Salmon's output via 
#' the files output by file_prep_Salmon.sh
#' 
#' @param cnt A raw count df with WORMBASE IDs as rownames and 16 columns:
#' 1 columns per replicate
#' @param geneIDs_in_both list of WORMBASE IDs present in both introns and exons raw counts
#' @param path_table path to a look up table wite gene names, WORMBASE IDs and chromosomes
#' @param chrom select a chromosome (ex:X) or all
#' @return raw counts with gene of interests named
#' @examples
#' cntEx <- get_names(cntEx, genes.in.both, gene_table, chrom="all")
get_names <- function(cnt, geneIDs_in_both, path_table, chrom="all") {
  
  cnt <- cnt[rownames(cnt) %in% geneIDs_in_both,]
  
  ### Import Gene name/chromosome look up table
  gene.table <- read.delim2(path_table)
  rownames(gene.table) <- gene.table[,1]
  
  if(chrom!="all") {
    gene.table <- gene.table[gene.table$Chromosome.scaffold.name==chrom,]
  }
  
  ### Replace gene IDs with gene names from the table
  ### Rows with Gene_IDs that aren't found are discarded
  ge <- gene.table$Gene.name[match(rownames(cnt), rownames(gene.table))]
  cnt <- cnt[!is.na(ge),]
  rownames(cnt) <- na.omit(ge)
  
  return(cnt)
}


#### Manipulations ####
#' From normalized counts, build a df with counts averaged among replicates
#' 
#' @param cnt A normalized count df with gene names as rownames and 16 columns:
#' 1 columns per replicate
#' @param logBefore Boolean, log2 transform the data before continuing
#' @return An eisa df with gene names as rownames and 6 columns:
#'  deltaExon, deltaIntron and their error bars 
#' @examples
#' get_deltas(cnt,logBefore=TRUE)
get_deltas <- function(cnt,logBefore=TRUE) {
  
  if(logBefore) {cnt <- log2(cnt)}
  
  cntEx.366 <- rowMeans(cnt[,1:4])
  cntEx.382 <- rowMeans(cnt[,5:8])
  cntIn.366 <- rowMeans(cnt[,9:12])
  cntIn.382 <- rowMeans(cnt[,13:16])
  
  varEx.366 <- RowVar(cnt[,1:4])
  varEx.382 <- RowVar(cnt[,5:8])
  varIn.366 <- RowVar(cnt[,9:12])
  varIn.382 <- RowVar(cnt[,13:16])
  
  covEx <- RowCov(cnt[,1:4],cnt[,5:8])
  covIn <- RowCov(cnt[,9:12],cnt[,13:16])
  
  sd.dIntron <- sqrt(((cntIn.382/cntIn.366)^2)*((varIn.366/(cntIn.366^2))+(varIn.382/(cntIn.382^2))-2*(covIn/(cntIn.366*cntIn.382))))
  sd.dExon <- sqrt(((cntEx.382/cntEx.366)^2)*((varEx.366/(cntEx.366^2))+(varEx.382/(cntEx.382^2))-2*(covEx/(cntEx.366*cntEx.382))))
  
  dExon <- cntEx.382/cntEx.366
  dIntron <- cntIn.382/cntIn.366
  
  if(logBefore) {
    sd.dIntron <- sqrt(varIn.366+varIn.382)
    sd.dExon <- sqrt(varEx.366+varEx.382)
    
    dExon <- cntEx.382-cntEx.366
    dIntron <- cntIn.382-cntIn.366
  }
  
  delta.cnt <- data.frame(dExon=dExon,
                          dIntron=dIntron,
                          lowI=dIntron-sd.dIntron/2,
                          highI=dIntron+sd.dIntron/2,
                          lowE=dExon-sd.dExon/2,
                          highE=dExon+sd.dExon/2)
  
  if(logBefore) {return(na.omit(delta.cnt[is.finite(rowSums(delta.cnt)),]))} 
  else {return(na.omit(log2(delta.cnt[is.finite(rowSums(log2(delta.cnt))),])))}
  
}

#' From normalized counts, build a df with counts averaged among replicates
#' 
#' @param cnt A normalized count df with gene names as rownames and 16 columns:
#' 1 columns per replicate
#' @return An averaged df with gene names as rownames and 12 columns:
#' cntEx366, cntIn366, cntEx382, cntIn382 and their error bars 
#' @examples
#' get_means(cnt)
get_means <- function(cnt) {
  
  cntEx.366 <- rowMeans(cnt[,1:4])
  cntEx.382 <- rowMeans(cnt[,5:8])
  cntIn.366 <- rowMeans(cnt[,9:12])
  cntIn.382 <- rowMeans(cnt[,13:16])
  
  sdEx.366 <- RowSD(cnt[,1:4])
  sdEx.382 <- RowSD(cnt[,5:8])
  sdIn.366 <- RowSD(cnt[,9:12])
  sdIn.382 <- RowSD(cnt[,13:16])
  
  delta.cnt <- data.frame(cntEx.366,cntEx.382,
                          cntIn.366,cntIn.382,
                          lowI366=cntIn.366-sdIn.366/2,
                          highI366=cntIn.366+sdIn.366/2,
                          lowI382=cntIn.382-sdIn.382/2,
                          highI382=cntIn.382+sdIn.382/2,
                          lowE366=cntEx.366-sdEx.366/2,
                          highE366=cntEx.366+sdEx.366/2,
                          lowE382=cntEx.382-sdEx.382/2,
                          highE382=cntEx.382+sdEx.382/2)
  return(na.omit(log2(delta.cnt[is.finite(rowSums(log2(delta.cnt))),])))
}
#### Plots ####

#' Build bar plot of DeltaIntron/Exon for each gene with error bars
#' 
#' @param delta.cnt An eisa df with gene names as rownames and 6 columns:
#'  deltaExon, deltaIntron and their error bars 
#' @return A bar plot of DeltaIntron/Exon for each gene (log2)
#' @examples
#' plot_col_deltas(delta.cnt)
plot_col_deltas <- function(delta.cnt) {
  
  delta.cnt <- cbind(rownames(delta.cnt),delta.cnt)
  colnames(delta.cnt)[1] <- "Genes"
  
  exoo <- dplyr::select(delta.cnt,-c(dIntron, lowI, highI)) %>% 
    pivot_longer(cols=c(dExon),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowE, high=highE)
  
  introo <- dplyr::select(delta.cnt,-c(dExon, lowE, highE)) %>% 
    pivot_longer(cols=c(dIntron),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowI, high=highI)
  
  eisa <- bind_rows(exoo, introo)
  
  eisa %>% 
    ggplot(aes(Genes, y = log2, fill = feature)) +
    geom_col(position = "dodge", width=.6, alpha=0.85) + theme_light() +
    geom_errorbar(aes(ymin = low, ymax = high),width=.6, size=0.2, position = "dodge") + 
    ggtitle('WT → dpy26cs\n(366 → 382)') +
    theme(plot.title = element_text(hjust=0.5, size=15, face="bold"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_fill_manual(name ="", values = c("dIntron" = "red4", "dExon" = "dodgerblue4"))
}

#' Build bar plot of counts of Intron/Exon for each gene with error bars
#' 
#' @param mean.cnt An averaged df with gene names as rownames and 12 columns created by get_means(): 
#' cntEx366, cntIn366, cntEx382, cntIn382 and their error bars 
#' @return A bar plot of counts for each gene (log2)
#' @examples
#' plot_col_means(mean.cnt)
plot_col_means <- function(mean.cnt) {
  mean.cnt <- cbind(rownames(mean.cnt),mean.cnt)
  colnames(mean.cnt)[1] <- "Genes"
  
  e366 <- dplyr::select(mean.cnt,c(Genes, cntEx.366, lowE366, highE366)) %>% 
    pivot_longer(cols=c(cntEx.366),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowE366, high=highE366)
  
  e382 <- dplyr::select(mean.cnt,c(Genes, cntEx.382, lowE382, highE382)) %>% 
    pivot_longer(cols=c(cntEx.382),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowE382, high=highE382)
  
  i382 <- dplyr::select(mean.cnt,c(Genes, cntIn.382, lowI382, highI382)) %>% 
    pivot_longer(cols=c(cntIn.382),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowI382, high=highI382)
  
  i366 <- dplyr::select(mean.cnt,c(Genes, cntIn.366, lowI366, highI366)) %>% 
    pivot_longer(cols=c(cntIn.366),
                 names_to = "feature", values_to = "log2") %>%
    dplyr::rename(low=lowI366, high=highI366)
  
  eisa <- bind_rows(e366, e382, i382, i366)
  
  eisa %>% 
    ggplot(aes(Genes, y = log2, fill = feature)) +
    geom_col(position = "dodge", width=.6, alpha=0.85) + theme_light() +
    geom_errorbar(aes(ymin = low, ymax = high),width=.6, size=0.2, position = "dodge") + 
    ggtitle('WT → dpy26cs\n(366 → 382)') +
    theme(plot.title = element_text(hjust=0.5, size=15, face="bold"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_fill_manual(name ="", values = c(cntEx.366 = "dodgerblue4", cntEx.382 = "dodgerblue4", cntIn.366 = "red4", cntIn.382 = "red4"))
}

#' Build scatter plot of DeltaIntron/Exon 
#' 
#' @param delta.cnt An eisa df with gene names as rownames and at least 2 columns:
#' $dIntron and $dExon
#' @param delta.cnt.signi A subset eisa df with gene names as rownames and at least 2 columns:
#' $dIntron and $dExon
#' @return A scatter plot of DeltaIntron/Exon with Pearson correlation coeff
#' @examples
#' scatter_deltas(delta.cnt,delta.cnt.signi)
scatter_deltas <- function(delta.cnt,delta.cnt.signi) {
  corEvI <- round(cor(delta.cnt$dIntron,delta.cnt$dExon), 3)
  corEvI.signi <- round(cor(delta.cnt.signi$dIntron,delta.cnt.signi$dExon), 3)
  
  `%notin%` <- Negate(`%in%`)
  delta.cnt <- delta.cnt[(rownames(delta.cnt) %notin% rownames(delta.cnt.signi)),]
  
  ggplot() + 
    geom_point(data=delta.cnt, mapping=aes(x=dIntron, y=dExon), alpha=0.5) +
    geom_point(data=delta.cnt.signi, mapping=aes(x=dIntron, y=dExon, color = "FDR<0.05\nfor dIntron or dExon"), alpha=0.75) +
    ggtitle(paste0('R = ',corEvI), subtitle = paste0('R = ',corEvI.signi)) +
    scale_colour_manual(name = NULL, values = c("FDR<0.05\nfor dIntron or dExon" = "red")) +
    theme_light() +
    theme(plot.title=element_text(size=12, face="italic", margin = margin(t=40, b = -38)),
          plot.subtitle=element_text(size=12, face="italic", color="red", margin = margin(t=40, b = -35)))# +lims(x = c(-2, 2), y = c(-2, 2))
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