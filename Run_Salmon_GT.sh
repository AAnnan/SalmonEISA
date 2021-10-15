#! /bin/bash

## Resource Allocation
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu
#SBATCH --mem=196G
#SBATCH –-cpus-per-task=24

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Salmon"

source ${CONDA_ACTIVATE} salmon

fa366b1_1=/scratch/aannan/rnaseq/data/366b1/366B1_L1_R1_001_F6TKUgD2sWF3.fastq.gz
fa366b1_2=/scratch/aannan/rnaseq/data/366b1/366B1_L2_R1_001_Aekf4XsT6VYC.fastq.gz

fa382b1_1=/scratch/aannan/rnaseq/data/382b1/382B1_L1_R1_001_bBV52jLO2VX2.fastq.gz
fa382b1_2=/scratch/aannan/rnaseq/data/382b1/382B1_L2_R1_001_wMojKBI3VV5t.fastq.gz

fa366b2_1=/scratch/aannan/rnaseq/data/366b2/366B2_L1_R1_001_VukPzfB7UznO.fastq.gz
fa366b2_1=/scratch/aannan/rnaseq/data/366b2/366B2_L2_R1_001_6eYI0PJPBw8Q.fastq.gz

fa366b3_1=/scratch/aannan/rnaseq/data/366b3/366B3_L1_R1_001_9YNWmzpY47EG.fastq.gz
fa366b3_2=/scratch/aannan/rnaseq/data/366b3/366B3_L2_R1_001_v1BuQK0d8xsQ.fastq.gz

fa366b4_1=/scratch/aannan/rnaseq/data/366b4/366B4_L1_R1_001_oOmUjBH4ogzp.fastq.gz
fa366b4_2=/scratch/aannan/rnaseq/data/366b4/366B4_L2_R1_001_KluTH7GqLaqT.fastq.gz

fa382b2_1=/scratch/aannan/rnaseq/data/382b2/382B2_L1_R1_001_tt73nYLLPlmX.fastq.gz
fa382b2_2=/scratch/aannan/rnaseq/data/382b2/382B2_L2_R1_001_u3AK2wvq6gvq.fastq.gz

fa382b3_1=/scratch/aannan/rnaseq/data/382b3/382B3_L1_R1_001_yMMdM9dtpDE9.fastq.gz
fa382b3_2=/scratch/aannan/rnaseq/data/382b3/382B3_L2_R1_001_JZWqrZz0eUf3.fastq.gz

fa382b4_1=/scratch/aannan/rnaseq/data/382b4/382B4_L1_R1_001_Y6Xklx9GlZkI.fastq.gz
fa382b4_2=/scratch/aannan/rnaseq/data/382b4/382B4_L2_R1_001_EGVBFP42rXm4.fastq.gz

cat c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa CE_genes_seq.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentrome.fa

salmon index -t CE_gentrome.fa -i EISA_index --decoys decoys.txt -k 19
#the k size selected here will act as the minimum acceptable length for a valid match.

salmon quant -i EISA_index -l A -r ${fa366b1_1} ${fa366b1_2} --validateMappings --writeMappings=366b1.sam --dumpEq -p 20 -o salmon_out/366b1 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${fa382b1_1} ${fa382b1_2} --validateMappings --writeMappings=382b1.sam --dumpEq -p 20 -o salmon_out/382b1 --seqBias --gcBias --numBootstraps 100 

salmon quant -i EISA_index -l A -r ${fa366b2_1} ${fa366b2_2} --validateMappings --writeMappings=366b2.sam --dumpEq -p 20 -o salmon_out/366b2 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${fa382b2_1} ${fa382b2_2} --validateMappings --writeMappings=382b2.sam --dumpEq -p 20 -o salmon_out/382b2 --seqBias --gcBias --numBootstraps 100  

salmon quant -i EISA_index -l A -r ${fa366b3_1} ${fa366b3_2} --validateMappings --writeMappings=366b3.sam --dumpEq -p 20 -o salmon_out/366b3 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${fa382b3_1} ${fa382b3_2} --validateMappings --writeMappings=382b3.sam --dumpEq -p 20 -o salmon_out/382b3 --seqBias --gcBias --numBootstraps 100  

salmon quant -i EISA_index -l A -r ${fa366b4_1} ${fa366b4_2} --validateMappings --writeMappings=366b4.sam --dumpEq -p 20 -o salmon_out/366b4 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${fa382b4_1} ${fa382b4_2} --validateMappings --writeMappings=382b4.sam --dumpEq -p 20 -o salmon_out/382b4 --seqBias --gcBias --numBootstraps 100  


conda deactivate
