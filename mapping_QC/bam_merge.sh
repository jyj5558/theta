#!/bin/bash
#SBATCH --job-name=bam_merge
#SBATCH -A fnrgenetics
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load bamtools

# cd $SLURM_SUBMIT_DIR

rm -rf autosomes.bam
rm -rf autosomes.sam
rm -rf out.bam
rm -rf out.sam

# merge bams
# list your SRA numbers in order to merge all bam files into one bam file
ls *.sam > samfiles.txt
sed "s/sam/bam/g" samfiles.txt > files.bamlist 

bamtools merge -list files.bamlist -out out.bam

# END