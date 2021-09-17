# SalmonEISA
Exon/Intron split Analysis (EISA) pipeline for RNA-Seq *C. elegans* data using GenomicRanges and Salmon

## Index

- [Prerequisites](https://github.com/AAnnan/SalmonEISA/#prerequisites)
- [Usage](https://github.com/AAnnan/SalmonEISA/#usage)
- [Output](https://github.com/AAnnan/SalmonEISA/#output)

## Prerequisites

A Salmon installation:
* Follow [this](https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon) or [this](https://salmon.readthedocs.io/en/latest/building.html).

The following R libraries: 
* `GenomicFeatures, tidyr, dplyr, edgeR, rtracklayer, BSgenome.Celegans.UCSC.ce11` 

The following files:
* *C. elegans* genome file in fasta format: **`c_elegans.PRJNA13758.WS279.genomic.fa`**
* *C. elegans* annotation file in GFF3 or GTF format: **`c_elegans.PRJNA13758.WS279.annotations.gff3`**
* *C. elegans* gene look-up table containing WORMBASE IDs, gene names and chromosomes. Find it [here](http://parasite.wormbase.org/biomart/martview?VIRTUALSCHEMANAME=parasite_mart&ATTRIBUTES=wbps_gene.default.feature_page.wbps_gene_id|wbps_gene.default.feature_page.external_gene_id|wbps_gene.default.feature_page.chromosome_name&FILTERS=wbps_gene.default.naive_filters.species_id_1010."caelegprjna13758"&VISIBLEPANEL=attributepanel). Click on `Results` on then `Go`.

## Usage

Run each script in this order:

1. **`Build_Exon_Intron_fastas.R`**
1. **`Run_Salmon.sh`**
1. **`Build_Rawcounts_EISA.sh`**
1. **`SalmonEISA.R`**

## Output

* Barplot showing **ΔExon** and **ΔIntron** for each gene (error bars represent the SD to the mean).
![Barplot](https://i.imgur.com/zixY90M.png)
* Scatterplot showing **ΔExon** in relation to **ΔIntron** for all genes and their Pearson correlation coefficient.
![Scatterplot](https://i.imgur.com/sQFQMTq.png)
