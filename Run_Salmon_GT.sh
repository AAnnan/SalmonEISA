#! /bin/bash

## Resource Allocation
#SBATCH --time=4-00:00:00
#SBATCH --partition=local
#SBATCH --mem=96G
#SBATCH –-cpus-per-task=20

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Salmon"

source ${CONDA_ACTIVATE} salmon

cat c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa CE_genes_seq.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentrome.fa

salmon index -t CE_gentrome.fa -i EISA_index --decoys decoys.txt -k 19
#the k size selected here will act as the minimum acceptable length for a valid match.

########################################################################################################################
########################################################################################################################
data=/scratch/aannan/Kruezi2013/data

#366 b 2020
faN2_GROseq_ctrl_1=${data}/GSM1056279_GRO-seq_N2_Emb/SRR639125.fastq
faN2_GROseq_ctrl_2=${data}/GSM1056279_GRO-seq_N2_Emb/SRR639126.fastq

faGROseq_sdc2_1=${data}/GSM1056281_GRO-seq_y93.sdc2RNAi_Emb/SRR639129.fastq
faGROseq_sdc2_2=${data}/GSM1056281_GRO-seq_y93.sdc2RNAi_Emb/SRR639130.fastq

########################################################################################################################
########################################################################################################################

#366 b 2020
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrl_1} --validateMappings -p 20 -o salmon_out/N2_GROseq_ctrl_1 --seqBias --gcBias --numBootstraps 100
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrl_2} --validateMappings -p 20 -o salmon_out/N2_GROseq_ctrl_2 --seqBias --gcBias --numBootstraps 100

salmon quant -i EISA_index -l A -r ${faGROseq_sdc2_1} --validateMappings -p 20 -o salmon_out/GROseq_sdc2_1 --seqBias --gcBias --numBootstraps 100
salmon quant -i EISA_index -l A -r ${faGROseq_sdc2_2} --validateMappings -p 20 -o salmon_out/GROseq_sdc2_2 --seqBias --gcBias --numBootstraps 100

conda deactivate




