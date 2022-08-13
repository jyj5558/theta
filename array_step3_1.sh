#!/bin/bash
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

#############################################
# this script will create:
# reference genome dictionary
# reference genome index files
# cleaned_sralist for the downstream analysis
#############################################

########################
#DO NOT EDIT BELOW CODE#
########################

# define variable
genus_species=$1

# load modules
module load bioinfo
module load bwa
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bamtools

# make a new directory to house processed alignment files
mkdir -p /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/

# move to reference directory
cd /scratch/bell/dewoody/theta/${genus_species}/*_ref/

# create dictionary for realignment
PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

# move to cleaned fastq files in preparation for alignment
cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/

# index reference 
bwa index -a bwtsw ../../*_ref/ref.fa

# capture all cleaned fastq files with variable
ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq > cleaned_sralist

# END
