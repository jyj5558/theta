#!/bin/bash
#SBATCH --job-name=S2_Genus-species
#SBATCH -A fnrdewoody
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=50G
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="your email address (e.g.abc123@purdue.edu) without quotation marks"

module load biocontainers
module load trim-galore
module load cutadapt
module load fastqc
module load sra-tools

####notes and usage####
#
##Notes##
#Edit the step2_SRA_download_clean.sh script, to add target species "genus-species" and the number of cpus "n" you allocated in the SBATCH command above
#script was written to be run with SLURM job scheduler
#
#Example:
#genus_species=Marmota-marmota-marmota
#n=32
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/abc/xyz/)
#usage:
#sbatch /scratch/bell/$USER/theta/step2_SRA_download_clean.sh
#
#Genus-species: this is used in the directory naming as Erangi suggested, to make browsing 
#a bit more doable for us humans
#
####end usage and notes####

################################
#Enter undefined variable below#
################################

genus_species=
n=

#########################
#DO NOT EDIT BELOW CODE #
#########################

#Extract releveant info from metadata file in theta directory for target species and save in species folder:

cd /path/to/theta/${genus_species}/
cat $CLUSTER_SCRATCH/theta/SRA_metadata/${genus_species}.txt | sed 's/ /_/g'  > ${genus_species}_SRA.txt


####create directories and download SRAs####
mkdir -p ./sra/raw
mkdir ./sra/cleaned
mkdir ./sra/aligned

cd ./sra/raw/

cat ../../${genus_species}_SRA.txt | cut -f 1 | tail -n +2 | while read g
do
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

# merge output from fastqc and check for FAILs
cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo "raw"
echo "$FAIL"
rm -r ${g}_?_fastqc* 

#trim adapters
trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

# check quality of trimmed reads
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
