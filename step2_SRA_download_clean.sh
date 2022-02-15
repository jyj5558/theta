#!/bin/bash
#SBATCH --job-name=TrimGalore
#SBATCH -A fnrdewoody
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=50G

module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc
module load sra-toolkit

####usage and notes####
#
#usage:
#step2_SRA_download_clean.sh Genus-species 
#Genus-species: this is used in the directory naming as Erangi suggested, to make browsing 
#a bit more doable for us humans
#
#Note: SRAS.txt must be in /scratch/bell/dewoody/theta/${genus_species}/
#SRAs.txt: comma separated list of SRAs to download (should be on a single row)
#
####end usage and notes####

genus_species=Marmota-marmota-marmota

#Extract releveant info from metadata file in theta directory for target species and save in species folder:
cd /scratch/bell/dewoody/theta/${genus_species}/
cat $CLUSTER_SCRATCH/theta/SRA_metadata/${genus_species}.txt | sed 's/ /_/g'  > ${genus_species}_SRA.txt


####create directories and download SRAs####
mkdir -p ./sra/raw
mkdir ./sra/cleaned
mkdir ./sra/aligned

cd /scratch/bell/dewoody/theta/${genus_species}/sra/raw/

cat ../../${genus_species}_SRA.txt | cut -f 1 | tail -n +2 | while read g
do
mkdir ${g}
cd ${g}
echo ${g} | xargs prefetch --max-size 500GB ${g} -O ./
echo ${g}.sra | xargs fasterq-dump -e 6
rm -f ${g}.sra

# check quality of reads
fastqc ${g}_1.fastq --extract --quiet
fastqc ${g}_2.fastq --extract --quiet
rm ${g}_1_fastqc.zip
rm ${g}_2_fastqc.zip

merge output from fastqc and check for FAILs
cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
FILE=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo mv ${g}_1_fastqc*  ${g}
echo mv ${g}_2_fastqc*  ${g}

#trim adapters
trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

# check quality of trimmed reads
cd ../cleaned
cat ${g}.sra_1.fastq_trimming_report.txt ${g}.sra_2.fastq_trimming_report.txt > ${g}_fastqc_summary.txt
rm ${g}.sra_1_val_1_fastqc.zip
rm ${g}.sra_2_val_2_fastqc.zip
FILE=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo "cleaned"
echo "$FILE"
rm ${g}_fastqc_summary.txt
rm ${g}_1_val_1_fastqc.html
rm ${g}_2_val_2_fastqc.html
done

# END
