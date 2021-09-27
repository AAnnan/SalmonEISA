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

####Extract the non-overlapping exonic or intronic parts from a GTF 
####with exonicParts and intronicParts from R package "GenomicFeatures"

#linked.to.single.gene.only
#excludes (TRUE) or includes (FALSE) exons (or introns) not linked to a gene or linked to more than one gene.
CE_exons <- exonicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)
CE_introns <- intronicParts(CE_gtf_txdb, linked.to.single.gene.only=FALSE)

#Remove introns that overlap with any exon from any transcript isoform
#Also remove exons that overlap with any intron
CE_introns_NO <- subsetByOverlaps(CE_introns, CE_exons, invert = TRUE, ignore.strand=TRUE)
CE_exons_NO <- subsetByOverlaps(CE_exons, CE_introns, invert = TRUE, ignore.strand=TRUE)

#Verify that no overlap is left
#findOverlaps(CE_exons_NO, CE_introns_NO)

# Trimming operations
e=10L # length of the trim on either side of the feature
k <- 2 # minimum length of a feature (will be filtered out by salmon index later anyway)
# INFO: Percentage of features that are below the threshold k
paste0("introns below ",k,' bp: ',round((length(which(CE_introns_NO@ranges@width<=k))/length(CE_introns_NO))*100,1),"% ",length(which(CE_introns_NO@ranges@width<=k)),"/",length(CE_introns_NO))
paste0("exons below ",k,' bp: ',round((length(which(CE_exons_NO@ranges@width<=k))/length(CE_exons_NO))*100,1),"% ",length(which(CE_exons_NO@ranges@width<=k)),"/",length(CE_exons_NO))

### INTRON TRIMMING
# Trim intron regions by e=ten basepairs on both sides 
# to ensure that intronic reads close to the exon junctions 
# are not counted as exonic reads. EISA 2015
CE_introns_NO@ranges@start <- CE_introns_NO@ranges@start + e
CE_introns_NO@ranges@width <- CE_introns_NO@ranges@width - e

### EXON TRIMMING
# Trim exon regions by e=ten basepairs on both sides 
# to ensure that intronic reads close to the exon junctions 
# are not counted as exonic reads. EISA 2015
#CE_exons_NO@ranges@start <- CE_exons_NO@ranges@start + e
#CE_exons_NO@ranges@width <- CE_exons_NO@ranges@width - e

### EXON EXPANDING (bad)
# expand exon regions by e=ten basepairs on both sides
#CE_exons_NO@ranges@start <- CE_exons_NO@ranges@start - e
#CE_exons_NO@ranges@width <- CE_exons_NO@ranges@width + e
# bad solution
# It obviously adds plenty of overlap between exonic and intronic parts...

# Filter small features
CE_introns_NO <- CE_introns_NO[-which(CE_introns_NO@ranges@width<k),]
CE_exons_NO <- CE_exons_NO[-which(CE_exons_NO@ranges@width<k),]

#Overlaid distributions of exon/intron length
ggplot() + 
  geom_histogram(data = data.frame(exon_sizes = log2(CE_exons_NO@ranges@width)), aes(x = exon_sizes, fill = "Exons"), binwidth=0.05, alpha = 0.5) +
  geom_histogram(data = data.frame(intron_sizes = log2(CE_introns_NO@ranges@width)), aes(x = intron_sizes, fill = "Introns"), binwidth=0.05, alpha = 0.5) +
  scale_fill_manual(name ="", values = c("Exons" = "red", "Introns" = "blue")) + xlab("Length in log2(bp)")

#### Make the Multi-fastas
#Import ce11 genome
CE_genome <- BSgenome.Celegans.UCSC.ce11
#ensembl style
seqlevelsStyle(CE_genome) <- "Ensembl"

#Retrieve sequence data
CE_exons_seq <- getSeq(CE_genome,CE_exons_NO)
CE_introns_seq <- getSeq(CE_genome,CE_introns_NO)

#Build a short description line for the multifasta
#Just exon number (artificial) and associated gene ID
get_des <- function(int_ex_GR, int_or_ex) {
  idx <- seq(from = 1, to = length(int_ex_GR), by = 1)
  gene_ids <- unlist(int_ex_GR@elementMetadata$gene_id)
  
  Des <- paste0(int_or_ex, idx,"\t", gene_ids)
  
  return(Des)
}

#Add the descritpion line to the DNAseq object
CE_exons_seq@ranges@NAMES <- get_des(CE_exons_NO, int_or_ex="exon")
CE_introns_seq@ranges@NAMES <- get_des(CE_introns_NO, int_or_ex="intron")

#export multifastas
writeXStringSet(CE_exons_seq, "CE_exons_seq.fa", format="fasta", append=FALSE, compress=FALSE)
writeXStringSet(CE_introns_seq, "CE_introns_seq.fa", format="fasta", append=FALSE, compress=FALSE)

###NOT USED

#Build a long description line for the multifasta
get_description <- function(int_ex_GR, int_or_ex) {
  transc_ids <- as.list(int_ex_GR@elementMetadata$tx_name)
  transc_ids <- sapply(transc_ids, paste, collapse = ':')
  chrs <- as.list(int_ex_GR@seqnames)
  start_pos <- int_ex_GR@ranges@start
  end_pos <- int_ex_GR@ranges@start+int_ex_GR@ranges@width-1
  stran <- as.list(int_ex_GR@strand)
  stran <- replace(stran, stran=="+","1")
  stran <- replace(stran, stran=="-","-1")
  gene_ids <- int_ex_GR@elementMetadata$gene_id
  
  Des <- paste(int_or_ex, " chromosome:","WBcel235:", chrs,':',
               start_pos,':',end_pos,":", stran, 
               " gene:", gene_ids," transcript:",transc_ids,  sep = "")
  return(Des)
}
