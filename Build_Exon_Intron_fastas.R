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
CE_cds <- cds(CE_gtf_txdb, columns=c("gene_id","tx_id", "tx_name"), filter=NULL)
CE_introns <- intronicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)

#Remove introns that overlap with any cds 
#Also remove cds that overlap with any intron
CE_introns_NO <- subsetByOverlaps(CE_introns, CE_cds, invert = TRUE, ignore.strand=FALSE)
CE_cds_NO <- subsetByOverlaps(CE_cds, CE_introns, invert = TRUE, ignore.strand=FALSE)
#Verify that no overlap is left
#findOverlaps(CE_introns_NO, CE_cds_NO)

# Trimming operations
e=95L # length of the trim/expansion on either side of the feature

### INTRON EXPANDING 
# Trim intron regions by e=ten basepairs on both sides 
# to ensure that intronic reads close to the exon junctions 
# are not counted as exonic reads. EISA 2015
CE_introns_NO@ranges@start <- CE_introns_NO@ranges@start - e
CE_introns_NO@ranges@width <- CE_introns_NO@ranges@width + e

#Overlaid distributions of exon/intron length
ggplot() + 
  geom_histogram(data = data.frame(cds_sizes = log2(CE_cds_NO@ranges@width)), aes(x = cds_sizes, fill = "CDS"), binwidth=0.05, alpha = 0.5) +
  geom_histogram(data = data.frame(intron_sizes = log2(CE_introns_NO@ranges@width)), aes(x = intron_sizes, fill = "Introns"), binwidth=0.05, alpha = 0.5) +
  scale_fill_manual(name ="", values = c("CDS" = "red", "Introns" = "blue")) + xlab("Length in log2(bp)")

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
#ensembl style
seqlevelsStyle(CE_genome) <- "Ensembl"

#Retrieve sequence data
CE_cds_seq <- getSeq(CE_genome,CE_cds_NO)
CE_introns_seq <- getSeq(CE_genome,CE_introns_NO)

#Add the descritpion line (gene name) to the DNAseq object
CE_introns_seq@ranges@NAMES <- unlist(CE_introns_NO@elementMetadata$gene_id)
CE_cds_seq@ranges@NAMES <- unlist(CE_cds_NO@elementMetadata$gene_id)

#export multifastas
writeXStringSet(CE_cds_seq, "CE_cds_seq.fa", format="fasta", append=FALSE, compress=FALSE)
writeXStringSet(CE_introns_seq, "CE_introns_seq.fa", format="fasta", append=FALSE, compress=FALSE)
