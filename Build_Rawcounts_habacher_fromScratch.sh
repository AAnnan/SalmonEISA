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
d="salmon_out"
sample=("N2_mock_" "N2_rege1RNAi_")
rep=("1" "2")

### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

# Extract the gene/transcript number and its corresponding gene ID
grep '>' c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa | sed 's/>//g' | sed 's/gene=//g' > gene_tx_list.txt
grep '>' CE_genes_seq.fa | sed 's/>//g' > gene_gen_list.txt


### Build read count columns by sample / replicate
for s in ${sample[@]}; do
	for r in ${rep[@]}; do
		echo "Extracting read count from "${s}${r}

		echo ${s}${r} > NumReads_gen_${s}${r}.sf && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_gen_list.txt ${d}/${s}${r}/quant.sf >> NumReads_gen_${s}${r}.sf
		echo ${s}${r} > NumReads_tx_${s}${r}.sf && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}${r}/quant.sf >> NumReads_tx_${s}${r}.sf
	done
done

### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

echo "Gene_IDs" > order_WBGenes_GE && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/${s[0]}${r[0]}/quant.sf gene_gen_list.txt >> order_WBGenes_GE
echo "Gene_IDs" > order_WBGenes_Tx && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/${s[0]}${r[0]}/quant.sf gene_tx_list.txt >> order_WBGenes_Tx

# Paste together all the columns
paste order_WBGenes_Tx NumReads_tx_N2_mock_1.sf NumReads_tx_N2_mock_2.sf NumReads_tx_N2_rege1RNAi_1.sf NumReads_tx_N2_rege1RNAi_2.sf > rawcounts_transcript.txt
paste order_WBGenes_GE NumReads_gen_N2_mock_1.sf NumReads_gen_N2_mock_2.sf NumReads_gen_N2_rege1RNAi_1.sf NumReads_gen_N2_rege1RNAi_2.sf > rawcounts_gene.txt

# Clean up
rm g*_t*_list.txt g*_g*_list.txt
rm NumReads*.sf
rm order*
