#!/bin/bash
#SBATCH --job-name=TrimGalore
#SBATCH -A fnrquail
#SBATCH -t 300:00:00
#SBATCH -N 1
#SBATCH -n 1

#load modules

module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc


for i in $(ls -1 *fastq | sed 's/_[1,2].fastq/ /g' | uniq)
do
echo "checking quality scores of paired-fastq sample $i"

# check quality of reads
fastqc ${i}_1.fastq --extract --quiet
fastqc ${i}_2.fastq --extract --quiet

echo fastqc is done, now merging fastqc output for $i

# merge output from fastqc and check for FAILs
cat ${i}_1_fastqc/summary.txt ${i}_2_fastqc/summary.txt > ${i}_fastqc_summary.txt

FILE=$(grep "FAIL" ${i}_fastqc_summary.txt)
echo "$FILE"

trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${i}_1.fastq ${i}_2.fastq

done


# remove excess files
rm -rf fastq_files.txt
rm -rf *_fastqc_summary.txt
rm -rf fail-txt

# END
