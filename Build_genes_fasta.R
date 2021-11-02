getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

MM_gtf_txdb <- makeTxDbFromGFF("mm_ensembl/Mus_musculus.GRCm38.102.gff3",format="gff3")

MM_genes <- genes(MM_gtf_txdb, columns="gene_id", filter=NULL, single.strand.genes.only=TRUE)

#### Make the Multi-fastas
#Import ce11 genome
MM_genome <- BSgenome.Mmusculus.UCSC.mm10
seqlevelsStyle(MM_genome) <- "ensembl" #ensembl style

#Retrieve sequence data
MM_genes_seq <- getSeq(MM_genome,MM_genes)

#Add the descritpion line to the DNAseq object
MM_genes_seq@ranges@NAMES <- MM_genes@elementMetadata$gene_id

#export Genes multifasta
writeXStringSet(MM_genes_seq, "MM_genes_seq.fa", format="fasta", append=FALSE, compress=FALSE)
