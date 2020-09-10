#!/bin/bash
#SBATCH --job-name=qcRef
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load bioawk
module load RepeatMasker/4.0.7
module load cmake/3.9.4

# cd $SLURM_SUBMIT_DIR

# make sure you have genmap installed and in your path
# https://github.com/cpockrandt/genmap
export PATH=/home/abruenic/genmap-build/bin:$PATH

# identify repeats in ref
RepeatMasker -species mammal ref.fa

# build an index of the fasta file(s) whose mappability you want to compute
genmap index -F ref.fa -I index -S 30

# compute mappability
# k = kmer of 100bp
# E = # two mismatches
genmap map -K 100 -E 2 -I index -O mappability -t -w -bg


# END