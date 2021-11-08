getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Celegans.UCSC.ce11)

# Make a TxDb GFF3 WS annotation
CE_gtf_txdb <- makeTxDbFromGFF("wormbase/c_elegans.PRJNA13758.WS279.annotations.gff3",format="gff3")

# Extract unspliced transcripts (=genes)
CE_genes <- genes(CE_gtf_txdb, columns="gene_id", filter=NULL, single.strand.genes.only=TRUE)

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
seqlevelsStyle(CE_genome) <- "Ensembl" #ensembl style

#Retrieve sequence data
CE_genes_seq <- getSeq(CE_genome,CE_genes)

#Add the descritpion line to the DNAseq object
CE_genes_seq@ranges@NAMES <- gsub("Gene:", "", CE_genes_seq@ranges@NAMES)

#export Genes multifasta
writeXStringSet(CE_genes_seq, "CE_genes_seq.fa", format="fasta", append=FALSE, compress=FALSE)
