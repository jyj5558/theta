#!/bin/bash
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

#####################################
# this script will create:
# sralist for the downstream analysis
# directories to store sra files
# script was written to be run with SLURM job scheduler
#####################################

########################
#DO NOT EDIT BELOW CODE#
########################

# define variables
genus_species=$1
THETA=/scratch/bell/dewoody/theta

# create directories and download SRAs
mkdir -p ./sra/raw
mkdir -p ./sra/cleaned
mkdir -p ./sra/aligned

# extract releveant info from metadata file in theta directory for target species and save in species folder:
cd ${THETA}/${genus_species}/
cat $CLUSTER_SCRATCH/theta/SRA_metadata/${genus_species}.txt | sed 's/ /_/g'  > ${genus_species}_SRA.txt
cat ${genus_species}_SRA.txt | cut -f 1 | tail -n +2 > ./sra/raw/sralist

# END
