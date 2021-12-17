#! /bin/bash

## Resource Allocation
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu
#SBATCH --mem=196G
#SBATCH –-cpus-per-task=20

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="STAR"

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

/home/aannan/soft/gffread/gffread /mnt/external.data/aannan/EISA/c_elegans.PRJNA13758.WS279.annotations.gff3  -T -o c_elegans.PRJNA13758.WS279.annotations.gtf

STAR --runMode genomeGenerate \
--genomeDir ./sequence --genomeFastaFiles c_elegans.PRJNA13758.WS279.genomic.fa \
--sjdbGTFfile c_elegans.PRJNA13758.WS279.annotations.gtf --runThreadN 20 \
--genomeSAindexNbases 12 --sjdbOverhang 49 

###########################################################################
###########################################################################

STAR --genomeDir ./sequence \
 --readFilesIn ${faAMA_1_1} ${faAMA_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/AMA_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faAMA_2_1} ${faAMA_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/AMA_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

###########################################################################
###########################################################################

STAR --genomeDir ./sequence \
 --readFilesIn ${faFLAVO_1_1} ${faFLAVO_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/FLAVO_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faFLAVO_2_1} ${faFLAVO_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/FLAVO_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

###########################################################################
###########################################################################

STAR --genomeDir ./sequence \
 --readFilesIn ${faN2_1_1} ${faN2_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/N2_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faN2_2_1} ${faN2_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/N2_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

###########################################################################
###########################################################################

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1m_1_1} ${faTIR1m_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1m_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1m_2_1} ${faTIR1m_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1m_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1m_3_1} ${faTIR1m_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1m_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1p_1_1} ${faTIR1p_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1p_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1p_2_1} ${faTIR1p_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1p_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTIR1p_3_1} ${faTIR1p_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TIR1p_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

###########################################################################
###########################################################################

 STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1m_1_1} ${faTOP1m_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1m_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1m_2_1} ${faTOP1m_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1m_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1m_3_1} ${faTOP1m_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1m_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1p_1_1} ${faTOP1p_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1p_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1p_2_1} ${faTOP1p_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1p_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP1p_3_1} ${faTOP1p_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP1p_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

###########################################################################
###########################################################################

 STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2m_1_1} ${faTOP2m_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2m_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2m_2_1} ${faTOP2m_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2m_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2m_3_1} ${faTOP2m_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2m_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2p_1_1} ${faTOP2p_1_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2p_1_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2p_2_1} ${faTOP2p_2_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2p_2_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

STAR --genomeDir ./sequence \
 --readFilesIn ${faTOP2p_3_1} ${faTOP2p_3_2} --readFilesCommand zcat \
 --outFileNamePrefix ./bamSTAR/TOP2p_3_ \
 --runThreadN 20 --outSAMmultNmax 1 \
 --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --alignIntronMax 5000 --quantMode GeneCounts --outWigType wiggle --outWigNorm RPM

conda deactivate

#sbatch --cpus-per-task 20 STAR

