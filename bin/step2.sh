#!/bin/bash

module purge    
module load biocontainers
module load trim-galore
module load cutadapt
module load fastqc
module load sra-tools

####notes and usage####
#
# This script downloads whole-genome resequencing data (SRA) of a focal population,
# and quality-control the data
# script was written to be run with SLURM job scheduler
#
####end of notes and usage####

########################
#DO NOT EDIT BELOW CODE#
########################

genus_species=$1
accession=$2
pathway=$3
assembly=$4
    
#Extract relevant info from metadata file in theta directory for target species and save in species folder:

cd $CLUSTER_SCRATCH/theta/${genus_species}/
cat $CLUSTER_SCRATCH/theta/SRA_metadata/${genus_species}.txt | sed 's/ /_/g' > ${genus_species}_SRA.txt


####create directories and download SRAs####
mkdir -p ./sra/raw
mkdir -p ./sra/cleaned
mkdir -p ./sra/aligned

cd ./sra/raw/

cat ../../${genus_species}_SRA.txt | cut -f 1 | tail -n +2 | while read g; do
  mkdir ${g}
  cd ${g}
  prefetch --max-size 500GB -O ./ ${g}
  fasterq-dump -e 32 --progress ${g}
  find . -name '*.fastq' -exec mv {} ../ \;
  cd ../
  rm -r ${g}

  fastqc ${g}_1.fastq --extract --quiet
  fastqc ${g}_2.fastq --extract --quiet
  rm ${g}_1_fastqc.zip
  rm ${g}_2_fastqc.zip

#Merge output from fastqc and check for FAILs
  cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
  FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
  echo "raw"
  echo "$FAIL"
  rm -r ${g}_?_fastqc* 

#Trim adapters
  trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

#Check quality of trimmed reads
  cd ../cleaned
  cat ${g}_1.fastq_trimming_report.txt ${g}_2.fastq_trimming_report.txt > ${g}_fastqc_summary.txt
  rm ${g}_1_val_1_fastqc.zip
  rm ${g}_2_val_2_fastqc.zip
  FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
  echo "cleaned"
  echo "$FAIL"
  rm ${g}_fastqc_summary.txt
  rm ${g}_1_val_1_fastqc.html
  rm ${g}_2_val_2_fastqc.html
  cd ../raw
done

#END                                                                                                                          