getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

#CE_gtf_txdb <- makeTxDbFromGFF("c_elegans.PRJNA13758.WS281.annotations.gff3",format="gff3")
CE_gtf_txdb <- makeTxDbFromGFF("wormbase/c_elegans.PRJNA13758.WS279.annotations.gff3",format="gff3")

CE_genes <- genes(CE_gtf_txdb, columns="gene_id", filter=NULL, single.strand.genes.only=TRUE)
CE_exons <- exonicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)
CE_introns <- intronicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)

#Remove introns that overlap with any exon from any transcript isoform
#Also remove exons that overlap with any intron
CE_introns_NO <- subsetByOverlaps(CE_introns, CE_exons, invert = TRUE, ignore.strand=TRUE)
CE_exons_NO <- subsetByOverlaps(CE_exons, CE_introns, invert = TRUE, ignore.strand=TRUE)

e=10L 
### EXON EXPANDING
# expand exon regions by e=ten basepairs on both sides
CE_exons_NO@ranges@start <- CE_exons_NO@ranges@start - e
CE_exons_NO@ranges@width <- CE_exons_NO@ranges@width + e

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
seqlevelsStyle(CE_genome) <- "Ensembl" #ensembl style

#Retrieve sequence data
CE_exons_seq <- getSeq(CE_genome,CE_exons_NO)
CE_genes_seq <- getSeq(CE_genome,CE_genes)

#Build a short description line for the multifasta
#Just exon number (artificial) and associated gene ID
get_des <- function(int_ex_GR, int_or_ex) {
  idx <- seq(from = 1, to = length(int_ex_GR), by = 1)
  gene_ids <- unlist(int_ex_GR@elementMetadata$gene_id)
  
  Des <- paste0(int_or_ex, idx,"\t", gene_ids)
  
  return(Des)
}
#Add the descritpion line to the DNAseq object
CE_genes_seq@ranges@NAMES <- CE_genes@elementMetadata$gene_id
CE_exons_seq@ranges@NAMES <- get_des(CE_exons_NO, int_or_ex="exon")

#export Genes multifasta
writeXStringSet(CE_exons_seq, "CE_exons_seq.fa", format="fasta", append=FALSE, compress=FALSE)
writeXStringSet(CE_genes_seq, "CE_genes_seq.fa", format="fasta", append=FALSE, compress=FALSE)
