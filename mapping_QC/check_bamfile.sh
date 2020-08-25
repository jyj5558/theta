#!/bin/bash
#SBATCH --job-name=check_bam
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
#module load picard-tools/2.9.0
module load samtools


# cd $SLURM_SUBMIT_DIR

# check bam file
#PicardCommandLine ValidateSamFile I=out.bam MODE=SUMMARY

# count number of reads in bam file
#samtools flagstat out.bam

#samtools flagstat corrected.bam

#samtools flagstat sorted.bam

# sort bamfiles
#sed "s/.bam//g" files.bamlist > bamfiles

#cat bamfiles| while read -r LINE

#do

#samtools sort -o ${LINE}_sorted.bam ${LINE}.bam

#done

# make list of sorted bams
#ls *_sorted.bam > sorted.bamlist

#samtools merge -b sorted.bamlist merged.bam

#samtools sort -o merged_sorted.bam merged.bam

samtools flagstat realigned_reads.bam

# END