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
# Salmon input (multifastas) and output (salmon_out)

# Variable names for file paths
d="../salmon_out"
rep1=("1" "2")
rep2=("1" "2" "3")

rAMA=()
rFLAVO=()
rN2=()
for r in ${rep1[@]}; do
	rAMA+="AMA_${r} "
	rFLAVO+="FLAVO_${r} "
	rN2+="N2_${r} "
done


### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

mkdir tmp
cd tmp
# Extract the gene/transcript number and its corresponding gene ID
grep '>' ../c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa | sed 's/>//g' | sed 's/gene=//g' > gene_tx_list.txt
grep '>' ../CE_genes_seq.fa | sed 's/>//g' > gene_gen_list.txt

echo "Gene_IDs" > order_WBGenes_GE && awk 'NR==FNR{a[$1];next} $1 in a{print $1}' ${d}/AMA_1/quant.sf gene_gen_list.txt >> order_WBGenes_GE
echo "Gene_IDs" > order_WBGenes_Tx && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/AMA_1/quant.sf gene_tx_list.txt >> order_WBGenes_Tx

mkdir gen
mkdir tx

### Build read count columns by sample / replicate
for s in ${r366b[@]} ${r382b[@]} ; do
	echo "Extracting read count from "${s}

	echo ${s} > gen/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_gen_list.txt ${d}/${s}/quant.sf >> gen/${s}
	echo ${s} > tx/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}/quant.sf >> tx/${s}
done

### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

# Contrasts: 
#WT2020_vs_dpy26cs2020
#${r366b[@]} ${r382b[@]}

# Paste together all the columns
cd gen
paste ../order_WBGenes_GE ${r366b[@]} ${r382b[@]} > ../order_WT2020_vs_dpy26cs2020_rawcounts_gene.txt

cd ../tx
paste ../order_WBGenes_Tx ${r366b[@]} ${r382b[@]} > ../WT2020_vs_dpy26cs2020_rawcounts_transcript.tx

# Clean up
mv *raw* ..
cd ..
rm order*
rm -r tmp