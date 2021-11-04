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
# EISA_23_09_2021_GT with Salmon input (multifastas) and outpt (salmon_out)

# Variable names for file paths
d="../salmon_out"

rna_sample=(N2_RNAseq_ctrl_1 N2_RNAseq_ctrl_2 N2_RNAseq_ctrl_3 N2_RNAseq_ctrl_4 RNAseq_sdc2_1 RNAseq_sdc2_2 RNAseq_sdc2_3)
gro_sample=(N2_GROseq_ctrl_1 N2_GROseq_ctrl_2 N2_GROseq_ctrlRNAi_1 N2_GROseq_ctrlRNAi_2 GROseq_sdc2_1 GROseq_sdc2_2)

### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

mkdir tmp
cd tmp
# Extract the gene/transcript number and its corresponding gene ID
grep '>' ../c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa | sed 's/>//g' | sed 's/gene=//g' > gene_tx_list.txt
grep '>' ../CE_genes_seq.fa | sed 's/>//g' > gene_gen_list.txt

echo "Gene_IDs" > order_WBGenes_GE && awk 'NR==FNR{a[$1];next} $1 in a{print $1}' ${d}/${gro_sample[0]}/quant.sf gene_gen_list.txt >> order_WBGenes_GE
echo "Gene_IDs" > order_WBGenes_Tx && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/${gro_sample[0]}/quant.sf gene_tx_list.txt >> order_WBGenes_Tx

mkdir gen
mkdir tx

### Build read count columns by sample / replicate
for s in ${rna_sample[@]}; do
	echo "Extracting read count from "${s}

	echo ${s} > gen/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_gen_list.txt ${d}/${s}/quant.sf >> gen/${s}
	echo ${s} > tx/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}/quant.sf >> tx/${s}
done

for s in ${gro_sample[@]}; do
	echo "Extracting read count from "${s}

	echo ${s} > gen/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_gen_list.txt ${d}/${s}/quant.sf >> gen/${s}
	echo ${s} > tx/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}/quant.sf >> tx/${s}
done


### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

# Paste together all the columns
cd gen
paste ../order_WBGenes_GE ${rna_sample[@]} > ../order_rawcounts_gene.txt
paste ../order_WBGenes_GE ${gro_sample[@]} > ../order_rawcounts_gene_gro.txt

cd ../tx
paste ../order_WBGenes_Tx ${rna_sample[@]} > ../rawcounts_transcript.txt
paste ../order_WBGenes_Tx ${gro_sample[@]} > ../rawcounts_transcript_gro.txt

cd ..
sed 's/Gene://g' order_rawcounts_gene_gro.txt > rawcounts_gene_gro.txt
sed 's/Gene://g' order_rawcounts_gene.txt > rawcounts_gene.txt

mv raw* ..
# Clean up
cd ..
rm -r tmp
