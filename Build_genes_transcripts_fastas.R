getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")
options(scipen=999) #prevent sci notation

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(ggplot2)

#CE_gtf_txdb <- makeTxDbFromGFF("c_elegans.PRJNA13758.WS281.annotations.gff3",format="gff3")
CE_gtf_txdb <- makeTxDbFromGFF("wormbase/c_elegans.PRJNA13758.WS279.annotations.gff3",format="gff3")

CE_genes <- genes(CE_gtf_txdb, columns="gene_id", filter=NULL, single.strand.genes.only=TRUE)
CE_transcripts <- transcripts(CE_gtf_txdb, columns=c("gene_id","tx_id", "tx_name"), filter=NULL)

#Attempt to retrieve only the longest transcript per gene
#a <- CE_transcripts[order(width(CE_transcripts), decreasing = TRUE)]
#a <- a[!duplicated(a$gene_id@unlistData)]

# INFO: Percentage of features that are below the threshold k
k <- 50
paste0("genes below ",k,' bp: ',round((length(which(CE_genes@ranges@width<=k))/length(CE_genes))*100,1),"% ",length(which(CE_genes@ranges@width<=k)),"/",length(CE_genes))
paste0("transcripts below ",k,' bp: ',round((length(which(CE_transcripts@ranges@width<=k))/length(CE_transcripts))*100,1),"% ",length(which(CE_transcripts@ranges@width<=k)),"/",length(CE_transcripts))

#Overlaid distributions of genes and transcripts length
ggplot() + 
  geom_histogram(data = data.frame(genes_sizes = log2(CE_genes@ranges@width)), aes(x = genes_sizes, fill = "Genes"), binwidth=0.05, alpha = 0.3) +
  geom_histogram(data = data.frame(transcripts_sizes = log2(CE_transcripts@ranges@width)), aes(x = transcripts_sizes, fill = "Transcripts"), binwidth=0.05, alpha = 0.3) +
  scale_fill_manual(name ="", values = c("Transcripts" = "red", "Genes" = "blue")) + xlab("Length in log2(bp)")

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
#ensembl style
seqlevelsStyle(CE_genome) <- "Ensembl"

#Retrieve sequence data
CE_genes_seq <- getSeq(CE_genome,CE_genes)
CE_transcripts_seq <- getSeq(CE_genome,CE_transcripts)

#Build a short description line for the multifasta
#Just gene/transcript number (artificial) and associated gene ID
get_des <- function(g_t_GR, g_or_t) {
  idx <- seq(from = 1, to = length(g_t_GR), by = 1)
  gene_ids <- unlist(g_t_GR@elementMetadata$gene_id)
  
  Des <- paste0(g_or_t, idx,"\t", gene_ids)
  
  return(Des)
}

#Add the descritpion line to the DNAseq object
CE_genes_seq@ranges@NAMES <- get_des(CE_genes, g_or_t="gene")
CE_transcripts_seq@ranges@NAMES <- get_des(CE_transcripts, g_or_t="transcript")

#export multifastas
writeXStringSet(CE_genes_seq, "CE_genes_seq.fa", format="fasta", append=FALSE, compress=FALSE)
writeXStringSet(CE_transcripts_seq, "CE_transcripts_seq.fa", format="fasta", append=FALSE, compress=FALSE)
