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

rAMA=()
rFLAVO=()
rN2=()
for r in ${rep1[@]}; do
	rAMA+="AMA_${r} "
	rFLAVO+="FLAVO_${r} "
	rN2+="N2_${r} "
done

rep2=("1" "2" "3")

rTIR1m=()
rTOP1m=()
rTOP2m=()
rTIR1p=()
rTOP1p=()
rTOP2p=()

for r2 in ${rep2[@]}; do
	rTIR1m+="TIR1m_${r2} "
	rTOP1m+="TOP1m_${r2} "
	rTOP2m+="TOP2m_${r2} "
	rTIR1p+="TIR1p_${r2} "
	rTOP1p+="TOP1p_${r2} "
	rTOP2p+="TOP2p_${r2} "
done

### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

mkdir tmp
cd tmp
# Extract the gene/transcript number and its corresponding gene ID
grep '>' ../c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa | sed 's/>//g' | sed 's/gene=//g' > gene_tx_list.txt

echo "Gene_IDs" > order_WBGenes_Tx && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/AMA_1/quant.sf gene_tx_list.txt >> order_WBGenes_Tx

mkdir tx

### Build read count columns by sample / replicate
for s in ${rAMA[@]} ${rFLAVO[@]} ${rN2[@]} ${rTIR1m[@]} ${rTOP1m[@]} ${rTOP2m[@]} ${rTIR1p[@]} ${rTOP1p[@]} ${rTOP2p[@]} ; do
	echo "Extracting read count from "${s}

	echo ${s} > tx/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}/quant.sf >> tx/${s}
done

### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

# Contrasts: 
#N2_vs_AMA
#${rN2[@]} ${rAMA[@]}
#N2_vs_FLAVO
#${rN2[@]} ${rFLAVO[@]}
#N2_vs_TIR1m
#${rN2[@]} ${rTIR1m[@]}
#N2_vs_TIR1p
#${rN2[@]} ${rTIR1p[@]}
#TOP1m_vs_TOP1p
#${rTOP1m[@]} ${rTOP1p[@]}
#TOP2m_vs_TOP2p
#${rTOP2m[@]} ${rTOP2p[@]}
#TIR1m_vs_TIR1p
#${rTIR1m[@]} ${rTIR1p[@]}

# Paste together all the columns

cd tx
paste ../order_WBGenes_Tx ${rN2[@]} ${rAMA[@]} > ../N2_vs_AMA_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rN2[@]} ${rFLAVO[@]} > ../N2_vs_FLAVO_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rN2[@]} ${rTIR1m[@]} > ../N2_vs_TIR1m_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rN2[@]} ${rTIR1p[@]} > ../N2_vs_TIR1p_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rTOP1m[@]} ${rTOP1p[@]} > ../TOP1m_vs_TOP1p_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rTOP2m[@]} ${rTOP2p[@]} > ../TOP2m_vs_TOP2p_rawcounts_tx_RNASEQ.txt
paste ../order_WBGenes_Tx ${rTIR1m[@]} ${rTIR1p[@]} > ../TIR1m_vs_TIR1p_rawcounts_tx_RNASEQ.txt

# Clean up
cd ../..
mkdir rnaseq_rawc
mv tmp/*raw* rnaseq_rawc/.
rm -r tmp