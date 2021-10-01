# SalmonEISA
Exon/Intron Split Analysis (EISA) pipeline for RNA-Seq *C. elegans* data using GenomicRanges and Salmon

## Index

- [Prerequisites](https://github.com/AAnnan/SalmonEISA/#prerequisites)
- [Usage](https://github.com/AAnnan/SalmonEISA/#usage)
- [Output](https://github.com/AAnnan/SalmonEISA/#output)

## Prerequisites

A Salmon installation:
* Follow [this](https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon) or [this](https://salmon.readthedocs.io/en/latest/building.html).

The following R libraries: 
* `GenomicRanges, GenomicFeatures, tidyr, dplyr, edgeR, rtracklayer, BSgenome.Celegans.UCSC.ce11` 

The following files:
* Cleaned, preprocessed reads in fasta format.
* *C. elegans* genome file in fasta format: **`c_elegans.PRJNA13758.WS279.genomic.fa`**
* *C. elegans* annotation file in GFF3 or GTF format: **`c_elegans.PRJNA13758.WS279.annotations.gff3`**
* *C. elegans* gene look-up table containing WORMBASE IDs, gene names and chromosomes. Find it [here](http://parasite.wormbase.org/biomart/martview?VIRTUALSCHEMANAME=parasite_mart&ATTRIBUTES=wbps_gene.default.feature_page.wbps_gene_id|wbps_gene.default.feature_page.external_gene_id|wbps_gene.default.feature_page.chromosome_name&FILTERS=wbps_gene.default.naive_filters.species_id_1010."caelegprjna13758"&VISIBLEPANEL=attributepanel). Click on `Results` on then `Go`.

## Usage

Take a look at the headers of each of the following scripts, modify them according to your data and run in order:

1. **`Build_Exon_Intron_fastas.R`**
1. **`Run_Salmon.sh`**
1. **`Build_Rawcounts_EISA.sh`**
1. **`SalmonEISA.R`**

## Output
Barplot showing **ΔExon** and **ΔIntron** for each gene (error bars represent the SD to the mean)             |  Scatterplot showing **ΔExon** in relation to **ΔIntron** for all genes (+ the Pearson correlation coefficient)
:-------------------------------------:|:-------------------------------------:
![](https://i.imgur.com/CpaFVJX.png)  |  ![](https://i.imgur.com/qvHMSq7.png)
