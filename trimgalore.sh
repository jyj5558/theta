#!/bin/bash
#SBATCH --job-name=TrimGalore
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load TrimGalore
module load cutadapt/2.5
module load fastqc

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script is checking for the quality of the SRAs
# if needed it remove adaptors and low quality sequence from reads and 	
# the following files to use in downstream analyses are created:	  
# *1_val_1.fq.gz and *_2_val_2.fq.gz									
#########################################################################

# check if your files are zipped
ls *sra*.fastq > unzipped_fastq_files.txt

if [ -f "unzipped_fastq_files.txt" ]
then

sed -i "s/.sra_1.fastq//g" unzipped_fastq_files.txt
sed -i "s/.sra_2.fastq//g" unzipped_fastq_files.txt

# this makes a list with the unique SRA names
sort -u unzipped_fastq_files.txt > unzipped.txt

cat unzipped.txt | while read -r LINE

do

gzip ${LINE}.sra_1.fastq
mv ${LINE}_1.fastq.gz

gzip ${LINE}.sra_2.fastq
mv ${LINE}_2.fastq.gz

done

fi

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

# check quality of reads
fastqc ${LINE}_1.fastq.gz
fastqc ${LINE}_2.fastq.gz

unzip ${LINE}_1_fastqc.zip
unzip ${LINE}_2_fastqc.zip

# merge output from fastqc and check for FAILs
cat ${LINE}_1_fastqc/summary.txt ${LINE}_2_fastqc/summary.txt > ${LINE}_fastqc_summary.txt

FILE=$(grep "FAIL" ${LINE}_fastqc_summary.txt < fail.txt)

if [ -s "$FILE" ]
then
# -q Trim low-quality ends from reads in addition to adapter removal. phred 20
# run fastqc without binning --fastqc_args
trim_galore --stringency 1 --length 30 --quality 20 \
--fastqc_args "--nogroup" -o $LINE --paired $READ1 $READ2
else
 	echo "$FILE OK"
fi

done

# remove excess files
rm -rf fastq_files.txt
rm -rf *_fastqc_summary.txt
rm -rf fail-txt

# END
