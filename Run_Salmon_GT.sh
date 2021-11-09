#! /bin/bash

## Resource Allocation
#SBATCH --time=2-00:00:00
#SBATCH --partition=local
#SBATCH --mem=160G
#SBATCH --cpus-per-task=30

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Salmon"

source ${CONDA_ACTIVATE} salmon

cat CE_genes_seq.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentrome.fa

salmon index -p 30 -t CE_gentrome.fa -i EISA_index --decoys decoys.txt -k 19
#the k size selected here will act as the minimum acceptable length for a valid match.

########################################################################################################################
########################################################################################################################
data=/scratch/aannan/GROseq_RNAseq/data

#Starved L1 GROseq
faGROseq_StarvedL1_1=${data}/Kruezi2013/GSM1056282_GRO-seq_N2_StarvedL1/SRR639131.fastq
faGROseq_StarvedL1_2=${data}/Kruezi2013/GSM1056282_GRO-seq_N2_StarvedL1/SRR639132.fastq

#Starved L1 RNAseq
faRNAseq_StarvedL1_1=${data}/Maxwell2012/12hr_Starved/SRR353593.fastq
faRNAseq_StarvedL1_2=${data}/Maxwell2012/12hr_Starved/SRR353596.fastq
faRNAseq_StarvedL1_3=${data}/Maxwell2012/12hr_Starved/SRR353600.fastq

#CTRL 1 GROseq
faN2_GROseq_ctrl_1=${data}/Kruezi2013/GSM1056279_GRO-seq_N2_Emb/SRR639125.fastq
faN2_GROseq_ctrl_2=${data}/Kruezi2013/GSM1056279_GRO-seq_N2_Emb/SRR639126.fastq

#CTRL RNAi 2 GROseq
faN2_GROseq_ctrlRNAi_1=${data}/Kruezi2013/GSM1056280_GRO-seq_controlRNAi_Emb/SRR639127.fastq
faN2_GROseq_ctrlRNAi_2=${data}/Kruezi2013/GSM1056280_GRO-seq_controlRNAi_Emb/SRR639128.fastq

#SDC2_RNAi GROseq
faGROseq_sdc2_1=${data}/Kruezi2013/GSM1056281_GRO-seq_y93.sdc2RNAi_Emb/SRR639129.fastq
faGROseq_sdc2_2=${data}/Kruezi2013/GSM1056281_GRO-seq_y93.sdc2RNAi_Emb/SRR639130.fastq

#CTRL RNAseq
faN2_RNAseq_ctrl_1=${data}/Crane2015/GSE59715_RNAseq_N2_Emb/SRR1523642.fastq
faN2_RNAseq_ctrl_2=${data}/Crane2015/GSE59715_RNAseq_N2_Emb/SRR1523643.fastq
faN2_RNAseq_ctrl_3=${data}/Crane2015/GSE59715_RNAseq_N2_Emb/SRR1523644.fastq
faN2_RNAseq_ctrl_4=${data}/Crane2015/GSE59715_RNAseq_N2_Emb/SRR1523645.fastq

#SDC2_RNAi RNAseq
faRNAseq_sdc2_1=${data}/Crane2015/GSE59715_RNAseq_y93.sdc2RNAi_Emb/SRR1523646.fastq
faRNAseq_sdc2_2=${data}/Crane2015/GSE59715_RNAseq_y93.sdc2RNAi_Emb/SRR1523647.fastq
faRNAseq_sdc2_3=${data}/Crane2015/GSE59715_RNAseq_y93.sdc2RNAi_Emb/SRR1523648.fastq

########################################################################################################################
########################################################################################################################

#Starved L1 GROseq
#salmon quant -i EISA_index -l A -r ${faGROseq_StarvedL1_1} --validateMappings -p 30 -o salmon_out/GROseq_StarvedL1_1 --seqBias --gcBias
#salmon quant -i EISA_index -l A -r ${faGROseq_StarvedL1_2} --validateMappings -p 30 -o salmon_out/GROseq_StarvedL1_2 --seqBias --gcBias

#Starved L1 RNAseq
#salmon quant -i EISA_index -l A -r ${faRNAseq_StarvedL1_1} --validateMappings -p 30 -o salmon_out/RNAseq_StarvedL1_1 --seqBias --gcBias
#salmon quant -i EISA_index -l A -r ${faRNAseq_StarvedL1_2} --validateMappings -p 30 -o salmon_out/RNAseq_StarvedL1_2 --seqBias --gcBias
#salmon quant -i EISA_index -l A -r ${faRNAseq_StarvedL1_3} --validateMappings -p 30 -o salmon_out/RNAseq_StarvedL1_3 --seqBias --gcBias

#CTRL 1 GROseq
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrl_1} --validateMappings -p 30 -o salmon_out/N2_GROseq_ctrl_1 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrl_2} --validateMappings -p 30 -o salmon_out/N2_GROseq_ctrl_2 --seqBias --gcBias

#CTRL RNAi 2 GROseq
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrlRNAi_1} --validateMappings -p 30 -o salmon_out/N2_GROseq_ctrlRNAi_1 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faN2_GROseq_ctrlRNAi_2} --validateMappings -p 30 -o salmon_out/N2_GROseq_ctrlRNAi_2 --seqBias --gcBias

#SDC2_RNAi GROseq
salmon quant -i EISA_index -l A -r ${faGROseq_sdc2_1} --validateMappings -p 30 -o salmon_out/GROseq_sdc2_1 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faGROseq_sdc2_2} --validateMappings -p 30 -o salmon_out/GROseq_sdc2_2 --seqBias --gcBias

#CTRL RNAseq
salmon quant -i EISA_index -l A -r ${faN2_RNAseq_ctrl_1} --validateMappings -p 30 -o salmon_out/N2_RNAseq_ctrl_1 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faN2_RNAseq_ctrl_2} --validateMappings -p 30 -o salmon_out/N2_RNAseq_ctrl_2 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faN2_RNAseq_ctrl_3} --validateMappings -p 30 -o salmon_out/N2_RNAseq_ctrl_3 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faN2_RNAseq_ctrl_4} --validateMappings -p 30 -o salmon_out/N2_RNAseq_ctrl_4 --seqBias --gcBias

#SDC2_RNAi RNAseq
salmon quant -i EISA_index -l A -r ${faRNAseq_sdc2_1} --validateMappings -p 30 -o salmon_out/RNAseq_sdc2_1 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faRNAseq_sdc2_2} --validateMappings -p 30 -o salmon_out/RNAseq_sdc2_2 --seqBias --gcBias
salmon quant -i EISA_index -l A -r ${faRNAseq_sdc2_3} --validateMappings -p 30 -o salmon_out/RNAseq_sdc2_3 --seqBias --gcBias

conda deactivate




