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

N2_mock_1=/home/aannan/Habacher2016/fastq/N2_mock_1/407F1-1_141014_D00404_0104_AC5HHAACXX_TGACCA-NoIndex_L006_R1_001.fastq.gz
N2_mock_2=/home/aannan/Habacher2016/fastq/N2_mock_2/407F2-1_141014_D00404_0104_AC5HHAACXX_ACAGTG-NoIndex_L006_R1_001.fastq.gz

N2_rege1RNAi_1=/home/aannan/Habacher2016/fastq/N2_rege1RNAi_1/407F3-1_141014_D00404_0104_AC5HHAACXX_TAGCTT-NoIndex_L006_R1_001.fastq.gz
N2_rege1RNAi_2=/home/aannan/Habacher2016/fastq/N2_rege1RNAi_2/407F4-1_141014_D00404_0104_AC5HHAACXX_GGCTAC-NoIndex_L006_R1_001.fastq.gz


cat c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa CE_genes_seq.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentrome.fa

#salmon index -t CE_gentrome.fa -i EISA_index --decoys decoys.txt -k 19
#the k size selected here will act as the minimum acceptable length for a valid match.

salmon quant -i EISA_index -l A -r ${N2_rege1RNAi_1} --validateMappings --writeMappings=N2_rege1RNAi_1.sam --dumpEq -p 20 -o salmon_out/N2_rege1RNAi_1 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${N2_rege1RNAi_2} --validateMappings --writeMappings=N2_rege1RNAi_2.sam --dumpEq -p 20 -o salmon_out/N2_rege1RNAi_2 --seqBias --gcBias --numBootstraps 100 

salmon quant -i EISA_index -l A -r ${N2_mock_1} --validateMappings --writeMappings=N2_mock_1.sam --dumpEq -p 20 -o salmon_out/N2_mock_1 --seqBias --gcBias --numBootstraps 100  
salmon quant -i EISA_index -l A -r ${N2_mock_2} --validateMappings --writeMappings=N2_mock_2.sam --dumpEq -p 20 -o salmon_out/N2_mock_2 --seqBias --gcBias --numBootstraps 100  


conda deactivate
