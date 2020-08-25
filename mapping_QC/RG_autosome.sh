#!/bin/bash
#SBATCH --job-name=realign_auto
#SBATCH -A fnrgenetics
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1

# cd $SLURM_SUBMIT_DIR

#GenomeAnalysisTK -T IndelRealigner -R ref.fa -I marked.bam \
#-targetIntervals forIndelRealigner.intervals -o realigned_reads.bam

# change RG in bam
#PicardCommandLine AddOrReplaceReadGroups \
#      I=realigned_reads.bam \
#      O=RG.bam \
#      RGID=group1 \
#      RGLB=lib1 \
#      RGPL=illumina \
#      RGPU=unit1 \
#      RGSM=sample1
#      "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1"

# remove the last line from the bed file. 
head -n -1 autosomes.bed > temp.txt ; mv temp.txt autosomes.bed

# https://broadinstitute.github.io/picard/explain-flags.html
samtools view -bh -f 3 -L autosomes.bed RG.bam > autosomes.bam

# sort bam file
samtools sort -T tmp autosomes.bam -o autosomes_sorted.bam

# convert to sam
samtools view -h -o autosomes.sam autosomes.bam

# check that the new bam is not corrupted
PicardCommandLine ValidateSamFile I=autosomes.sam IGNORE_WARNINGS=true \
MODE=SUMMARY O=autosomes_out.txt TMP_DIR=`pwd`/tmp

rm -rf autosomes_${LINE}.sam
rm -rf autosomes_${LINE}_sorted*

# END
