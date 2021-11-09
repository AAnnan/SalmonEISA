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

repb1=("1" "2" "3" "4")
repb2=("1" "2" "4" "5")
repCD=("1" "2" "3")

r366b=()
r382b=()
r775b=()

r784b=()

r366C=()
r828C=()
r844C=()
r366D=()
r821D=()
r822A=()
r822C=()
r823D=()

for r in ${repb1[@]}; do
	r366b+="366b${r} "
	r382b+="382b${r} "
	r775b+="775b${r} "
done

for r in ${repb2[@]}; do
	r784b+="784b${r} "
done

for r in ${repCD[@]}; do
	r366C+="366C${r} "
	r366D+="366D${r} "
	r821D+="821D${r} "
	r822A+="822A${r} "
	r822C+="822C${r} "
	r823D+="823D${r} "
	r828C+="828C${r}"
	r844C+="844C${r}"
done
### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

mkdir tmp
cd tmp
# Extract the gene/transcript number and its corresponding gene ID
grep '>' ../c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa | sed 's/>//g' | sed 's/gene=//g' > gene_tx_list.txt
grep '>' ../CE_genes_seq.fa | sed 's/>//g' > gene_gen_list.txt

echo "Gene_IDs" > order_WBGenes_GE && awk 'NR==FNR{a[$1];next} $1 in a{print $1}' ${d}/366C1/quant.sf gene_gen_list.txt >> order_WBGenes_GE
echo "Gene_IDs" > order_WBGenes_Tx && awk 'NR==FNR{a[$1];next} $1 in a{print $2}' ${d}/366C1/quant.sf gene_tx_list.txt >> order_WBGenes_Tx

mkdir gen
mkdir tx

### Build read count columns by sample / replicate
for s in ${r366b[@]} ${r382b[@]} ${r775b[@]} ${r784b[@]} ${r366C[@]} ${r366D[@]} ${r821D[@]} ${r822A[@]} ${r822C[@]} ${r828C[@]} ${r844C[@]} ${r823D[@]}; do
	echo "Extracting read count from "${s}

	echo ${s} > gen/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_gen_list.txt ${d}/${s}/quant.sf >> gen/${s}
	echo ${s} > tx/${s} && awk 'NR==FNR{a[$1];next} $1 in a{print $5}' gene_tx_list.txt ${d}/${s}/quant.sf >> tx/${s}
done

### Build final raw geneic/transcriptic counts
echo "Building raw geneic/transcriptic count files"

# Contrasts: 
#WT2020_vs_dpy26cs2020
#${r366b[@]} ${r382b[@]}
#WT2021_vs_dpy26cs2020
#${r366C[@]} ${r382b[@]}
#WT2021aWT2020_vs_dpy26cs2020
#${r366b[@]} ${r366C[@]} ${r382b[@]}
#WT2020_vs_kle2_2020
#${r366b[@]} ${r775b[@]}
#WT2020_vs_scc1_2020
#${r366b[@]} ${r784b[@]}

#TIR1wtA_vs_TIR1sd3degA
#${r823D[@]} ${r822A[@]}
#TIR1sd3deg_vs_TIR1sd3degA
#${r822C[@]} ${r822A[@]}
#TIR1wtA_vs_dpy26TIR1sd3degA
#${r823D[@]} ${r821D[@]}
#TIR1sd3deg_vs_dpy26TIR1sd3degA
#${r822C[@]} ${r821D[@]}
#WT2021_vs_coh1_2021
#${r366C[@]} ${r828C[@]}
#WT2021_vs_scc1coh1_2021
#${r366C[@]} ${r844C[@]}

# Paste together all the columns
cd gen
paste ../order_WBGenes_GE ${r366b[@]} ${r382b[@]} > ../WT2020_vs_dpy26cs2020_rawcounts_gene.txt
paste ../order_WBGenes_GE ${r366b[@]} ${r775b[@]} > ../WT2020_vs_kle2_2020_rawcounts_gene.txt
paste ../order_WBGenes_GE ${r366b[@]} ${r784b[@]} > ../WT2020_vs_scc1_2020_rawcounts_gene.txt
paste ../order_WBGenes_GE ${r366C[@]} ${r828C[@]} > ../WT2021_vs_coh1_2021_rawcounts_gene.txt
paste ../order_WBGenes_GE ${r366C[@]} ${r844C[@]} > ../WT2021_vs_scc1coh1_2021_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r366C[@]} ${r382b[@]} > ../WT2021_vs_dpy26cs2020_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r366b[@]} ${r366C[@]} ${r382b[@]} > ../WT2021aWT2020_vs_dpy26cs2020_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r823D[@]} ${r822A[@]} > ../TIR1wtA_vs_TIR1sd3degA_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r822C[@]} ${r822A[@]} > ../TIR1sd3deg_vs_TIR1sd3degA_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r823D[@]} ${r821D[@]} > ../TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_gene.txt
#paste ../order_WBGenes_GE ${r822C[@]} ${r821D[@]} > ../TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_gene.txt


cd ../tx
paste ../order_WBGenes_Tx ${r366b[@]} ${r382b[@]} > ../WT2020_vs_dpy26cs2020_rawcounts_transcript.txt
paste ../order_WBGenes_Tx ${r366b[@]} ${r775b[@]} > ../WT2020_vs_kle2_2020_rawcounts_transcript.txt
paste ../order_WBGenes_Tx ${r366b[@]} ${r784b[@]} > ../WT2020_vs_scc1_2020_rawcounts_transcript.txt
paste ../order_WBGenes_Tx ${r366C[@]} ${r828C[@]} > ../WT2021_vs_coh1_2021_rawcounts_transcript.txt
paste ../order_WBGenes_Tx ${r366C[@]} ${r844C[@]} > ../WT2021_vs_scc1coh1_2021_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r366C[@]} ${r382b[@]} > ../WT2021_vs_dpy26cs2020_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r366b[@]} ${r366C[@]} ${r382b[@]} > ../WT2021aWT2020_vs_dpy26cs2020_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r823D[@]} ${r822A[@]} > ../TIR1wtA_vs_TIR1sd3degA_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r822C[@]} ${r822A[@]} > ../TIR1sd3deg_vs_TIR1sd3degA_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r823D[@]} ${r821D[@]} > ../TIR1wtA_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt
#paste ../order_WBGenes_Tx ${r822C[@]} ${r821D[@]} > ../TIR1sd3deg_vs_dpy26TIR1sd3degA_rawcounts_transcript.txt


# Clean up
cd ../..
mv tmp/*rawcounts* .
rm -r tmp
