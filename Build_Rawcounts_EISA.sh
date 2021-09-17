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
# CE_introns_seq.fa, CE_exons_seq.fa, salmon_out

# Variable names for file paths
d="./salmon_out"
sample=("366" "382")
rep=("b1" "b2" "b3" "b4")

### Build Feature-Gene correspondance table
echo "Building Feature->GeneID correspondance table"

# Extract the exon/intron number and its corresponding gene ID
grep '>intron[0-9]\+' CE_introns_seq.fa | sed 's/>intron//g' > gene_int_list.txt
grep '>exon[0-9]\+' CE_exons_seq.fa | sed 's/>exon//g' > gene_ex_list.txt

# Sort the table by ascending exon/intron number
sort -k1 -o gene_int_list.txt gene_int_list.txt
sort -k1 -o gene_ex_list.txt gene_ex_list.txt

### Build read count columns by sample / replicate
for s in ${sample[@]}; do
	for r in ${rep[@]}; do
		echo "Extracting read count from "${s}${r}

		# Extract read counts from Salmon's output for introns and exons
		grep 'intron[0-9]\+' ${d}/${s}${r}/quant.sf | sed 's/intron//g' > quant_int_${s}${r}.sf
		grep 'exon[0-9]\+' ${d}/${s}${r}/quant.sf | sed 's/exon//g' > quant_ex_${s}${r}.sf

		# Sort by ascending exon/intron number
		sort -k1 -o quant_int_${s}${r}.sf quant_int_${s}${r}.sf
		sort -k1 -o quant_ex_${s}${r}.sf quant_ex_${s}${r}.sf

		# Select NumReads column (#5)
		cut -f5 quant_int_${s}${r}.sf > numReads_int_${s}${r}.sf
		cut -f5 quant_ex_${s}${r}.sf > numReads_ex_${s}${r}.sf
		
		# Add header
		cat <(echo ${s}${r}) numReads_int_${s}${r}.sf > int
		cat <(echo ${s}${r}) numReads_ex_${s}${r}.sf > ex
		cat int > numReads_int_${s}${r}.sf
		cat ex > numReads_ex_${s}${r}.sf
		#Clean up
		rm int ex
	done
done

### Build final raw exonic/intronic counts
echo "Building raw exonic/intronic count files"

# Since not all exons/introns present in the multifasta could be used (length, etc)
# This step filters out unused features (they're the same for each salmon run)
join -j 1 -t $'\t' -o 1.2 gene_int_list.txt quant_int_${s[0]}${r[0]}.sf > gen_int_list.txt
join -j 1 -t $'\t' -o 1.2 gene_ex_list.txt quant_ex_${s[0]}${r[0]}.sf > gen_ex_list.txt

# Add header
cat <(echo "Gene_IDs") gen_int_list.txt > ge_i_list.txt
cat <(echo "Gene_IDs") gen_ex_list.txt > ge_e_list.txt

# Paste together all the columns
paste ge_i_list.txt numReads_int_366b1.sf numReads_int_366b2.sf numReads_int_366b3.sf numReads_int_366b4.sf numReads_int_382b1.sf numReads_int_382b2.sf numReads_int_382b3.sf numReads_int_382b4.sf > rawcounts_intronic.txt
paste ge_e_list.txt numReads_ex_366b1.sf numReads_ex_366b2.sf numReads_ex_366b3.sf numReads_ex_366b4.sf numReads_ex_382b1.sf numReads_ex_382b2.sf numReads_ex_382b3.sf numReads_ex_382b4.sf > rawcounts_exonic.txt

# Clean up
rm g*_i*_list.txt g*_e*_list.txt
rm numReads*.sf
rm quant_*.sf

