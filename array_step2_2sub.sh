#!/bin/bash

#################################################################################################
# this script downloads sra files, remove low-quality reads, and trim adaptor sequences					#
# the following files to use in downstream analyses are created:	  	#
# *_1_val_1.fq and *_2_val_2.fq files for each sra file												#
#################################################################################################

#########################
#DO NOT EDIT BELOW CODE #
#########################

# load modules needed
module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc
module load sra-toolkit

# define variables
genus_species=$1
sra=$2
cpus_per_task=$3

echo "The task $SLURM_ARRAY_TASK_ID handles $sra"

# move to directory that will house files
cd /scratch/bell/dewoody/theta/${genus_species}/sra/raw

# remotely download sra fastq files 
mkdir ${sra}
cd ${sra}
echo ${sra} | xargs prefetch --max-size 500GB -O ./
echo ${sra}.sra | xargs fasterq-dump -e ${cpus_per_task} --progress
find . -name '*.fastq' -exec mv {} ../ \;
cd ../
rm -r ${sra}

# initial quality check of fastq files using fastqc
fastqc ${sra}_1.fastq --extract --quiet
fastqc ${sra}_2.fastq --extract --quiet
rm ${sra}_1_fastqc.zip
rm ${sra}_2_fastqc.zip

# merge output from fastqc and check for FAILs
cat ${sra}_1_fastqc/summary.txt ${sra}_2_fastqc/summary.txt > ${sra}_fastqc_summary.txt
FILE=$(grep "FAIL" ${sra}_fastqc_summary.txt)
echo "raw"
echo "$FILE"
rm -r ${sra}_?_fastqc* 

# trim adapters
trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${sra}_1.fastq ${sra}_2.fastq

# check quality of trimmed reads
cd ../cleaned
cat ${sra}_1.fastq_trimming_report.txt ${sra}_2.fastq_trimming_report.txt > ${sra}_fastqc_summary.txt
rm ${sra}_1_val_1_fastqc.zip
rm ${sra}_2_val_2_fastqc.zip
FILE=$(grep "FAIL" ${sra}_fastqc_summary.txt)
echo "cleaned"
echo "$FILE"
rm ${sra}_fastqc_summary.txt
rm ${sra}_1_val_1_fastqc.html
rm ${sra}_2_val_2_fastqc.html

# END
