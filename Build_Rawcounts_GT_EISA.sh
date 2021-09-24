#! /bin/bash

## Resource Allocation
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu
#SBATCH --mem=32G
#SBATCH –-cpus-per-task=4

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="file_prep_Salmon"

## Run in folder containing: 
# EISA_23_09_2021_Tx, EISA_23_09_2021_Genes with Salmon input and outpt

# Variable names for file paths
g="EISA_23_09_2021_Genes"
t="EISA_23_09_2021_Tx"
d="salmon_out"
sample=("366" "382")
rep=("b1" "b2" "b3" "b4")

### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

# Extract the gene/transcript number and its corresponding gene ID
grep '>transcript[0-9]\+' ${t}/CE_transcripts_seq.fa | sed 's/>transcript//g' > gene_tx_list.txt
grep '>gene[0-9]\+' ${g}/CE_genes_seq.fa | sed 's/>gene//g' > gene_gen_list.txt

# Sort the table by ascending gene/transcript number
sort -k1 -o gene_tx_list.txt gene_tx_list.txt
sort -k1 -o gene_gen_list.txt gene_gen_list.txt

### Build read count columns by sample / replicate
for s in ${sample[@]}; do
	for r in ${rep[@]}; do
		echo "Extracting read count from "${s}${r}

		# Extract read counts from Salmon's output for transcripts and genes
		grep 'transcript[0-9]\+' ${t}/${d}/${s}${r}/quant.sf | sed 's/transcript//g' > quant_tx_${s}${r}.sf
		grep 'gene[0-9]\+' ${g}/${d}/${s}${r}/quant.sf | sed 's/gene//g' > quant_gen_${s}${r}.sf

		# Sort by ascending gene/transcript number
		sort -k1 -o quant_tx_${s}${r}.sf quant_tx_${s}${r}.sf
		sort -k1 -o quant_gen_${s}${r}.sf quant_gen_${s}${r}.sf

		# Select numReads column (#5)
		cut -f5 quant_tx_${s}${r}.sf > NumReads_tx_${s}${r}.sf
		cut -f5 quant_gen_${s}${r}.sf > NumReads_gen_${s}${r}.sf
		
		# Add header
		cat <(echo ${s}${r}) NumReads_tx_${s}${r}.sf > tx
		cat <(echo ${s}${r}) NumReads_gen_${s}${r}.sf > gen
		cat tx > NumReads_tx_${s}${r}.sf
		cat gen > NumReads_gen_${s}${r}.sf
		#Clean up
		rm tx gen
	done
done

### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

# Since not all genes/transcripts present in the multifasta could be used by Salmon and output (length, etc)
# This step filters out unused features (they're the same for each salmon run)
join -j 1 -t $'\t' -o 1.2 gene_tx_list.txt quant_tx_${s[0]}${r[0]}.sf > gen_tx_list.txt
join -j 1 -t $'\t' -o 1.2 gene_gen_list.txt quant_gen_${s[0]}${r[0]}.sf > gen_gen_list.txt

# Add header
cat <(echo "Gene_IDs") gen_tx_list.txt > ge_t_list.txt
cat <(echo "Gene_IDs") gen_gen_list.txt > ge_g_list.txt

# Paste together all the columns
paste ge_t_list.txt NumReads_tx_366b1.sf NumReads_tx_366b2.sf NumReads_tx_366b3.sf NumReads_tx_366b4.sf NumReads_tx_382b1.sf NumReads_tx_382b2.sf NumReads_tx_382b3.sf NumReads_tx_382b4.sf > rawcounts_transcript.txt
paste ge_g_list.txt NumReads_gen_366b1.sf NumReads_gen_366b2.sf NumReads_gen_366b3.sf NumReads_gen_366b4.sf NumReads_gen_382b1.sf NumReads_gen_382b2.sf NumReads_gen_382b3.sf NumReads_gen_382b4.sf > rawcounts_gene.txt

# Clean up
rm g*_t*_list.txt g*_g*_list.txt
rm NumReads*.sf
rm *quant_*.sf

