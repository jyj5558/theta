#!/bin/bash
#SBATCH --job-name=filterCheck
#SBATCH -A fnrquail
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load gsl
module load zlib
module load gcc
module load samtools
module load bioawk

# cd $SLURM_SUBMIT_DIR

export PATH=/home/abruenic/angsd/:$PATH
export PATH=/home/abruenic/angsd/misc/:$PATH

# index bam and ref
samtools index realigned_reads.bam
samtools faidx ref.fa

# check filtered sites
output="/home/leopard/users/patricia/errorRates/filtered"
bam="/home/leopard/users/krishang/bam_filtered/results/bam/3241.bam"
ref="/home/leopard/users/patricia/errorRates/filtered/3241_filt_consensus.fa.gz"
anc="/home/leopard/users/fonseca/REF2/GCF_001857705.1_PanPar1.0_genomic.s.fna"
bamlist="/home/leopard/users/patricia/errorRates/filtered/filteredBamlist.txt"
SITES="/home/leopard/users/genis/regions/autosomes.filter_map_depth_inbreed.1K.regions"

angsd -i $goodbam -doCounts 1 -doFasta 2 -out $output/3241_filt_consensus

samtools faidx $ref 

angsd -doAncError 1 -anc $anc -ref $ref -nThreads 20 -sites $SITES -bam $bamlist -out $output/errorRate_filt_leopard

Rscript /programs/albrecht/angsd/R/estError.R file="errorRate_filt_leopard.ancError"