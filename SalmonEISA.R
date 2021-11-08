getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

##### Files (counts, geneTable) and parameters (EdgeR)
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"
txFile <- "rawcounts/ALL_GT/WT2020_vs_dpy26cs2020_rawcounts_transcript.txt"
geFile <- "rawcounts/ALL_GT/WT2020_vs_dpy26cs2020_rawcounts_gene.txt"

conditions <- c("WT","WT","WT","WT","dpy26","dpy26","dpy26","dpy26")
#######

SalmonEISA(geFile, txFile, gene_table, conditions)

