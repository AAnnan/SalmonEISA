#! /bin/bash

source ${CONDA_ACTIVATE} py38

gsutil -u stergachis cp gs://sra-pub-src-5/SRR2924297/407F1-1_141014_D00404_0104_AC5HHAACXX_TGACCA-NoIndex_L006_R1_001.fastq.gz ./N2_mock_1/.
gsutil -u stergachis cp gs://sra-pub-src-6/SRR2924298/407F2-1_141014_D00404_0104_AC5HHAACXX_ACAGTG-NoIndex_L006_R1_001.fastq.gz ./N2_mock_2/.
gsutil -u stergachis cp gs://sra-pub-src-4/SRR2924299/407F3-1_141014_D00404_0104_AC5HHAACXX_TAGCTT-NoIndex_L006_R1_001.fastq.gz ./N2_rege1RNAi_1/.
gsutil -u stergachis cp gs://sra-pub-src-8/SRR2924300/407F4-1_141014_D00404_0104_AC5HHAACXX_GGCTAC-NoIndex_L006_R1_001.fastq.gz ./N2_rege1RNAi_2/.

conda deactivate