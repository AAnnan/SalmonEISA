getwd()
setwd("~/Documents/Project_IZB/biodata/rnaseq/")

library(edgeR)
library(ggplot2)
source("/Users/aa/Documents/GitHub/SalmonEISA/SalmonEISA_func.R")

# input files and parameters
gene_table <- "wormbase/c_elegans.PRJNA13758.WS279.TableGeneIDs.tsv"

SalmonEISA(gene_table=gene_table,geFile="rawcounts/ALL_GT/WT2020_vs_dpy26cs2020_rawcounts_gene.txt",txFile="rawcounts/ALL_GT/WT2020_vs_dpy26cs2020_rawcounts_transcript.txt",conditions=c("WT","WT","WT","WT","dpy26","dpy26","dpy26","dpy26"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/REST_ALL/WT2020_vs_dpy26cs2020_rawcounts_gene.txt",txFile="rawcounts/REST_ALL/WT2020_vs_dpy26cs2020_rawcounts_transcript.txt",conditions=c("WT","WT","WT","WT","dpy26","dpy26","dpy26","dpy26"))

SalmonEISA(gene_table=gene_table,geFile="rawcounts/REST_ALL/WT2020_vs_kle2_2020_rawcounts_gene.txt",txFile="rawcounts/REST_ALL/WT2020_vs_kle2_2020_rawcounts_transcript.txt",conditions=c("WT","WT","WT","WT","kle2","kle2","kle2","kle2"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/REST_ALL/WT2020_vs_scc1_2020_rawcounts_gene.txt",txFile="rawcounts/REST_ALL/WT2020_vs_scc1_2020_rawcounts_transcript.txt",conditions=c("WT","WT","WT","WT","scc1","scc1","scc1","scc1"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/REST_ALL/WT2021_vs_coh1_2021_rawcounts_gene.txt",txFile="rawcounts/REST_ALL/WT2021_vs_coh1_2021_rawcounts_transcript.txt",conditions=c("WT","WT","WT","coh1","coh1","coh1"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/REST_ALL/WT2021_vs_scc1coh1_2021_rawcounts_gene.txt",txFile="rawcounts/REST_ALL/WT2021_vs_scc1coh1_2021_rawcounts_transcript.txt",conditions=c("WT","WT","WT","scc1coh1","scc1coh1","scc1coh1"))

SalmonEISA(gene_table=gene_table,geFile="rawcounts/ALL_GT/TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",txFile="rawcounts/ALL_GT/TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",conditions=c("TIR1sd3deg","TIR1sd3deg","TIR1sd3deg","dpy26TIR1sd3degA","dpy26TIR1sd3degA","dpy26TIR1sd3degA"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/ALL_GT/TIR1sd3deg_vs_TIR1sd3degA_rawcounts_gene.txt",txFile="rawcounts/ALL_GT/TIR1sd3deg_vs_TIR1sd3degA_rawcounts_transcript.txt",conditions=c("TIR1sd3deg","TIR1sd3deg","TIR1sd3deg","TIR1sd3degA","TIR1sd3degA","TIR1sd3degA"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/ALL_GT/TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",txFile="rawcounts/ALL_GT/TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",conditions=c("TIR1wtA","TIR1wtA","TIR1wtA","dpy26TIR1sd3degA","dpy26TIR1sd3degA","dpy26TIR1sd3degA"))
SalmonEISA(gene_table=gene_table,geFile="rawcounts/ALL_GT/TIR1wtA_vs_TIR1sd3degA_rawcounts_gene.txt",txFile="rawcounts/ALL_GT/TIR1wtA_vs_TIR1sd3degA_rawcounts_transcript.txt",conditions=c("TIR1wtA","TIR1wtA","TIR1wtA","TIR1sd3degA","TIR1sd3degA","TIR1sd3degA"))

#SalmonEISA(gene_table,"WT2021aWT2020_vs_dpy26cs2020_rawcounts_gene.txt","WT2021aWT2020_vs_dpy26cs2020_rawcounts_transcript.txt",c("WT","WT","WT","WT","WT","WT","WT","exp","exp","exp","exp"),c(1,1,1,1,1,1,1,2,2,2,2),"WT2021aWT2020-dpy26cs2020")
#SalmonEISA(gene_table,"WT2021_vs_dpy26cs2020_rawcounts_gene.txt","WT2021_vs_dpy26cs2020_rawcounts_transcript.txt",c("WT","WT","WT","exp","exp","exp","exp"),c(1,1,1,2,2,2,2),"WT2021-dpy26cs2020")



geFiles <- list("WT2020_vs_dpy26cs2020_rawcounts_gene.txt",
                "TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",
                "TIR1sd3deg_vs_TIR1sd3degA_rawcounts_gene.txt",
                "TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_gene.txt",
                "TIR1wtA_vs_TIR1sd3degA_rawcounts_gene.txt",
                "WT2021aWT2020_vs_dpy26cs2020_rawcounts_gene.txt",
                "WT2021_vs_dpy26cs2020_rawcounts_gene.txt")

txFiles <- list("WT2020_vs_dpy26cs2020_rawcounts_transcript.txt",
                "TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1sd3deg_vs_TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt",
                "TIR1wtA_vs_TIR1sd3degA_rawcounts_transcript.txt",
                "WT2021aWT2020_vs_dpy26cs2020_rawcounts_transcript.txt",
                "WT2021_vs_dpy26cs2020_rawcounts_transcript.txt")

conditions <- list(c("WT","WT","WT","WT","exp","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp"),
                   c("WT","WT","WT","WT","WT","WT","WT","exp","exp","exp","exp"),
                   c("WT","WT","WT","exp","exp","exp","exp"))

group_s <- list(c(1,1,1,1,2,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,1,1,1,1,2,2,2,2),
                c(1,1,1,2,2,2,2))

conds <- list("WT2020-dpy26cs2020",
              "TIR1sd3deg-dpy26TIR1sd3degA",
              "TIR1sd3deg-TIR1sd3degA",
              "TIR1wtA-dpy26TIR1sd3degA",
              "TIR1wtA-TIR1sd3degA",
              "WT2021aWT2020-dpy26cs2020",
              "WT2021-dpy26cs2020")
