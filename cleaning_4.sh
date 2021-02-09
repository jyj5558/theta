#!/bin/bash
#SBATCH --job-name=TrimGalore_4
#SBATCH -A fnrquail
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 1

#load modules

module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc



echo Moving to directory that houses raw paired-end fastq files for TWO channel chemistry
cd ./four_ch/raw

for g in $(ls -1 *fastq | sed 's/_[1,2].fastq/ /g' | uniq)
do
echo "checking quality scores of paired-fastq sample $g"

# check quality of reads
fastqc ${g}_1.fastq --extract --quiet
fastqc ${g}_2.fastq --extract --quiet

echo fastqc is done, now merging fastqc output for $g

# merge output from fastqc and check for FAILs
cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${i}_fastqc_summary.txt

FILE=$(grep "FAIL" ${g}_fastqc_summary.txt)
echo "$FILE"

trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

done


# remove excess files
rm -rf fastq_files.txt
rm -rf *_fastqc_summary.txt
rm -rf fail-txt

# END
