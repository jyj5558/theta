#!/bin/bash
#SBATCH --job-name=TrimGalore
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load TrimGalore
module load cutadapt/2.5

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script is remove adaptors and low quality sequence from reads  	#
# the following files to use in downstream analyses are created:	  	#
# *1_val_1.fq.gz and *_2_val_2.fq.gz									#
#########################################################################

# this makes a list of your fastq.gz read files 
ls *.fastq.gz > fastq_files.txt

# you might have to change the /_*.fastq.gz/ parts to fit with your reads
sed -i "s/_1.fastq.gz//g" fastq_files.txt
sed -i "s/_2.fastq.gz//g" fastq_files.txt

# this makes a list with the unique SRA names
sort -u fastq_files.txt > fastq.txt

cat fastq.txt | while read -r LINE

# again here you have to make sure that the _1/_2 part fit with your reads
READ1=${LINE}_1.fastq.gz
READ2=${LINE}_2.fastq.gz

do
# -q Trim low-quality ends from reads in addition to adapter removal. phred 20
# run fastqc without binning --fastqc_args
trim_galore --stringency 1 --length 30 --quality 20 --paired \
--fastqc_args "--nogroup" -o $LINE $READ1 $READ2

done

# remove excess files
rm -rf fastq_files.txt

# END
