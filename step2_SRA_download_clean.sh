#!/bin/bash
#SBATCH --job-name=S2_Genus-species
#SBATCH -A fnrdewoody
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=50G
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc
module load sra-toolkit

####notes and usage####
#
##Notes##
#Edit the step2_SRA_download_clean.sh script, to add target species "genus-species"
#Example:
#genus_species=Marmota-marmota-marmota
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/)
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

#########################
#DO NOT EDIT BELOW CODE #
#########################

#Extract releveant info from metadata file in theta directory for target species and save in species folder:

cd /scratch/bell/dewoody/theta/${genus_species}/
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
echo ${g} | xargs prefetch --max-size 500GB -O ./
echo ${g}.sra | xargs fasterq-dump -e 20 --progress
find . -name '*.fastq' -exec mv {} ../ \;
cd ../
rm -r ${g}

fastqc ${g}_1.fastq --extract --quiet
fastqc ${g}_2.fastq --extract --quiet
rm ${g}_1_fastqc.zip
rm ${g}_2_fastqc.zip

# merge output from fastqc and check for FAILs
cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
FILE=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo "raw"
echo "$FILE"
rm  ${g}_?_fastqc* 

#trim adapters
trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

# check quality of trimmed reads
cd ../cleaned
cat ${g}_1.fastq_trimming_report.txt ${g}_2.fastq_trimming_report.txt > ${g}_fastqc_summary.txt
rm ${g}_1_val_1_fastqc.zip
rm ${g}_2_val_2_fastqc.zip
FILE=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo "cleaned"
echo "$FILE"
rm ${g}_fastqc_summary.txt
rm ${g}_1_val_1_fastqc.html
rm ${g}_2_val_2_fastqc.html
cd ../raw
done
