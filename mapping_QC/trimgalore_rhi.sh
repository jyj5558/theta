#!/bin/bash
#SBATCH --job-name=TrimGalore
#SBATCH -A fnrdewoody
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools
module load TrimGalore

# cd $SLURM_SUBMIT_DIR

#sed -i "s/_1.fastq.gz//g" fastq_files.txt
#sed -i "s/_2.fastq.gz//g" fastq_files.txt

#sort -u fastq_files.txt > fastq.txt
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR943/SRR943144/SRR943144_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR943/SRR943144/SRR943144_2.fastq.gz

LINE=$"SRR943144"
READ1=$"$SRR943144_1.fastq.gz"
READ2=$"$SRR943144_2.fastq.gz"

# -q Trim low-quality ends from reads in addition to adapter removal. phred 20
# run fastqc without binning --fastqc_args
trim_galore --stringency 1 --length 30 --quality 20 --paired \
--fastqc_args "--nogroup" -o $LINE $READ1 $READ2


# EEND
