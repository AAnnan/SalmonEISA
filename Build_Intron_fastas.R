#this script will only create the intron list, get the transcript list from wormbase
getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")
options(scipen=999) #prevent sci notation

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(ggplot2)

CE_gtf_txdb <- makeTxDbFromGFF("wormbase/c_elegans.PRJNA13758.WS279.annotations.gff3",format="gff3")

#linked.to.single.gene.only
#excludes (TRUE) or includes (FALSE) exons (or introns) not linked to a gene or linked to more than one gene.
CE_exons <- exonicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)
CE_introns <- intronicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)

#Remove introns that overlap with any cds 
#Also remove cds that overlap with any intron
CE_introns_NO <- subsetByOverlaps(CE_introns, CE_exons, invert = TRUE, ignore.strand=FALSE)

# Trimming operations
e=99L # length of the trim/expansion on either side of the feature

### INTRON EXPANDING 
# Trim intron regions by e=ten basepairs on both sides 
# to ensure that intronic reads close to the exon junctions 
# are not counted as exonic reads. EISA 2015
CE_introns_NO@ranges@start <- CE_introns_NO@ranges@start - e
CE_introns_NO@ranges@width <- CE_introns_NO@ranges@width + e

#Overlaid distributions of exon/intron length
ggplot() + 
  geom_histogram(data = data.frame(exon_sizes = log2(CE_exons@ranges@width)), aes(x = exon_sizes, fill = "Exons"), binwidth=0.05, alpha = 0.5) +
  geom_histogram(data = data.frame(intron_sizes = log2(CE_introns_NO@ranges@width)), aes(x = intron_sizes, fill = "Introns"), binwidth=0.05, alpha = 0.5) +
  scale_fill_manual(name ="", values = c("Exons" = "red", "Introns" = "blue")) + xlab("Length in log2(bp)")

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
#ensembl style
seqlevelsStyle(CE_genome) <- "Ensembl"

#Retrieve sequence data
CE_introns_seq <- getSeq(CE_genome,CE_introns_NO)

#Build a short description line for the multifasta
#Just exon number (artificial) and associated gene ID
get_des <- function(int_ex_GR, int_or_ex) {
  idx <- seq(from = 1, to = length(int_ex_GR), by = 1)
  gene_ids <- unlist(int_ex_GR@elementMetadata$gene_id)
  
  Des <- paste0(int_or_ex, idx,"\t", gene_ids)
  
  return(Des)
}

#Can't use the gene names directly because Salmon doesn't take duplicates
#CE_introns_seq@ranges@NAMES <- unlist(CE_introns_NO@elementMetadata$gene_id)
CE_introns_seq@ranges@NAMES <- get_des(CE_introns_NO, int_or_ex="intron")

#export multifastas
writeXStringSet(CE_introns_seq, "CE_introns_seq.fa", format="fasta", append=FALSE, compress=FALSE)