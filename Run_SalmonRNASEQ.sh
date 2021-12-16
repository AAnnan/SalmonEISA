#! /bin/bash

## Resource Allocation
#SBATCH --time=3-00:00:00
#SBATCH --partition=gpu
#SBATCH --mem=340G
#SBATCH –-ntasks=72

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Salmon"

source ${CONDA_ACTIVATE} salmon
faAMA_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/AMA_1_1.fq.gz
faAMA_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/AMA_1_2.fq.gz
faAMA_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/AMA_2_1.fq.gz
faAMA_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/AMA_2_2.fq.gz


faFLAVO_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/FLAVO_1_1.fq.gz
faFLAVO_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/FLAVO_1_2.fq.gz
faFLAVO_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/FLAVO_2_1.fq.gz
faFLAVO_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/FLAVO_2_2.fq.gz


faN2_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/N2_1_1.fq.gz
faN2_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/N2_1_2.fq.gz
faN2_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/N2_2_1.fq.gz
faN2_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/N2_2_2.fq.gz

faTIR1m_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_1_1.fq.gz
faTIR1m_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_1_2.fq.gz
faTIR1m_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_2_1.fq.gz
faTIR1m_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_2_2.fq.gz
faTIR1m_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_3_1.fq.gz
faTIR1m_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_m_3_2.fq.gz


faTIR1p_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_1_1.fq.gz
faTIR1p_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_1_2.fq.gz
faTIR1p_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_2_1.fq.gz
faTIR1p_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_2_2.fq.gz
faTIR1p_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_3_1.fq.gz
faTIR1p_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TIR1_p_3_2.fq.gz


faTOP1m_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_1_1.fq.gz
faTOP1m_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_1_2.fq.gz
faTOP1m_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_2_1.fq.gz
faTOP1m_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_2_2.fq.gz
faTOP1m_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_3_1.fq.gz
faTOP1m_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_m_3_2.fq.gz


faTOP1p_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_1_1.fq.gz
faTOP1p_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_1_2.fq.gz
faTOP1p_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_2_1.fq.gz
faTOP1p_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_2_2.fq.gz
faTOP1p_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_3_1.fq.gz
faTOP1p_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP1_p_3_2.fq.gz


faTOP2m_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_1_1.fq.gz
faTOP2m_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_1_2.fq.gz
faTOP2m_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_2_1.fq.gz
faTOP2m_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_2_2.fq.gz
faTOP2m_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_3_1.fq.gz
faTOP2m_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_m_3_2.fq.gz


faTOP2p_1_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_1_1.fq.gz
faTOP2p_1_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_1_2.fq.gz
faTOP2p_2_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_2_1.fq.gz
faTOP2p_2_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_2_2.fq.gz
faTOP2p_3_1=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_3_1.fq.gz
faTOP2p_3_2=/scratch/aannan/16_12_eisa_bolaji/Bolaji20211216/TOP2_p_3_2.fq.gz


cat c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentromeRNASEQ.fa

salmon index -p 72 -t CE_gentromeRNASEQ.fa -i RNASEQ_index --decoys decoys.txt -k 17
#the k size selected here will act as the minimum acceptable length for a valid match.

#AMA
salmon quant -i RNASEQ_index -l A -1 ${faAMA_1_1} -2 ${faAMA_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/AMA_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faAMA_2_1} -2 ${faAMA_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/AMA_2 --seqBias --gcBias 

#FLAVO
salmon quant -i RNASEQ_index -l A -1 ${faFLAVO_1_1} -2 ${faFLAVO_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/FLAVO_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faFLAVO_2_1} -2 ${faFLAVO_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/FLAVO_2 --seqBias --gcBias 

#N2
salmon quant -i RNASEQ_index -l A -1 ${faN2_1_1} -2 ${faN2_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/N2_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faN2_2_1} -2 ${faN2_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/N2_2 --seqBias --gcBias 

#TIR1m
salmon quant -i RNASEQ_index -l A -1 ${faTIR1m_1_1} -2 ${faTIR1m_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1m_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTIR1m_2_1} -2 ${faTIR1m_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1m_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTIR1m_3_1} -2 ${faTIR1m_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1m_3 --seqBias --gcBias 

#TIR1p
salmon quant -i RNASEQ_index -l A -1 ${faTIR1p_1_1} -2 ${faTIR1p_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1p_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTIR1p_2_1} -2 ${faTIR1p_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1p_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTIR1p_3_1} -2 ${faTIR1p_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TIR1p_3 --seqBias --gcBias 

#TOP1m
salmon quant -i RNASEQ_index -l A -1 ${faTOP1m_1_1} -2 ${faTOP1m_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1m_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP1m_2_1} -2 ${faTOP1m_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1m_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP1m_3_1} -2 ${faTOP1m_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1m_3 --seqBias --gcBias 

#TOP1p
salmon quant -i RNASEQ_index -l A -1 ${faTOP1p_1_1} -2 ${faTOP1p_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1p_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP1p_2_1} -2 ${faTOP1p_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1p_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP1p_3_1} -2 ${faTOP1p_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP1p_3 --seqBias --gcBias 

#TOP2m
salmon quant -i RNASEQ_index -l A -1 ${faTOP2m_1_1} -2 ${faTOP2m_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2m_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP2m_2_1} -2 ${faTOP2m_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2m_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP2m_3_1} -2 ${faTOP2m_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2m_3 --seqBias --gcBias 

#TOP2p
salmon quant -i RNASEQ_index -l A -1 ${faTOP2p_1_1} -2 ${faTOP2p_1_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2p_1 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP2p_2_1} -2 ${faTOP2p_2_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2p_2 --seqBias --gcBias 
salmon quant -i RNASEQ_index -l A -1 ${faTOP2p_3_1} -2 ${faTOP2p_3_2} --validateMappings -p 72 -o salmon_out_RNASEQ/TOP2p_3 --seqBias --gcBias 

conda deactivate

#sbatch --cpus-per-task 72 SALM