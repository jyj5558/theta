#!/bin/bash
#SBATCH --job-name=SRAs
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

