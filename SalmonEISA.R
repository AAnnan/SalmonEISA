getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/N2_vs_AMA_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/N2_vs_AMA_rawcounts_transcript.txt",
           conditions=c("N2","N2","AMA","AMA"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/N2_vs_FLAVO_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/N2_vs_FLAVO_rawcounts_transcript.txt",
         conditions=c("N2","N2","FLAVO","FLAVO"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/N2_vs_TIR1m_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/N2_vs_TIR1m_rawcounts_transcript.txt",
           conditions=c("N2","N2","TIR1m","TIR1m","TIR1m"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/N2_vs_TIR1p_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/N2_vs_TIR1p_rawcounts_transcript.txt",
           conditions=c("N2","N2","TIR1p","TIR1p","TIR1p"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/TOP1m_vs_TOP1p_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/TOP1m_vs_TOP1p_rawcounts_transcript.txt",
           conditions=c("TOP1m","TOP1m","TOP1m","TOP1p","TOP1p","TOP1p"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/TOP2m_vs_TOP2p_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/TOP2m_vs_TOP2p_rawcounts_transcript.txt",
           conditions=c("TOP2m","TOP2m","TOP2m","TOP2p","TOP2p","TOP2p"))

EISA(gene_table=gene_table,geFile="rawcounts/bolaji1612/TIR1m_vs_TIR1p_rawcounts_gene.txt",
           txFile="rawcounts/bolaji1612/TIR1m_vs_TIR1p_rawcounts_transcript.txt",
           conditions=c("TIR1m","TIR1m","TIR1m","TIR1p","TIR1p","TIR1p"))
###########################################################################################
VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/N2_vs_AMA_rawcounts_tx_RNASEQ.txt",
     conditions=c("N2","N2","AMA","AMA"))

VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/N2_vs_FLAVO_rawcounts_tx_RNASEQ.txt",
     conditions=c("N2","N2","FLAVO","FLAVO"))

VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/N2_vs_TIR1m_rawcounts_tx_RNASEQ.txt",
     conditions=c("N2","N2","TIR1m","TIR1m","TIR1m"))

VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/N2_vs_TIR1p_rawcounts_tx_RNASEQ.txt",
     conditions=c("N2","N2","TIR1p","TIR1p","TIR1p"))

VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/TOP1m_vs_TOP1p_rawcounts_tx_RNASEQ.txt",
     conditions=c("TOP1m","TOP1m","TOP1m","TOP1p","TOP1p","TOP1p"))

VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/TOP2m_vs_TOP2p_rawcounts_tx_RNASEQ.txt",
     conditions=c("TOP2m","TOP2m","TOP2m","TOP2p","TOP2p","TOP2p"))
VOLC(gene_table=gene_table,
     txFile="rawcounts/bolaji1612_rnaseq/TIR1m_vs_TIR1p_rawcounts_tx_RNASEQ.txt",
     conditions=c("TIR1m","TIR1m","TIR1m","TIR1p","TIR1p","TIR1p"))
