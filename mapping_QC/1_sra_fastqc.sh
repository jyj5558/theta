#!/bin/sh

#PBS -N sra_QC
#PBS -q fnrdewoody
#PBS -l nodes=1:ppn=2,naccesspolicy=singleuser
#PBS -l walltime=300:00:00
#PBS -m abe

module purge
module load bioinfo
module load fastqc

cd $PBS_O_WORKDIR

# get sra
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR192/DRR192150/DRR192150_1.fastq.gz .
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR192/DRR192150/DRR192150_2.fastq.gz .

sra=$"DRR192150"

# this defines sample and SRA
#sra=$(echo "$LINE" | cut -d$'\t' -f2)
#echo $sra

FILE1=${sra}_1.fastq.gz
FILE2=${sra}_2.fastq.gz

# QC individual files
fastqc $FILE1
fastqc $FILE2

unzip ${sra}_1_fastqc.zip 
unzip ${sra}_2_fastqc.zip

grep "Per base sequence quality" ${sra}_1_fastqc/summary.txt >> temp_log.txt
grep "Adapter Content" ${sra}_1_fastqc/summary.txt >> temp_log.txt

grep "Per base sequence quality" ${sra}_2_fastqc/summary.txt >> temp_log.txt
grep "Adapter Content" ${sra}_2_fastqc/summary.txt >> temp_log.txt


if grep -q 'FAIL' temp_log.txt

then

# AdaptorRemoval
# trimming + output Phred 33 score
# https://github.com/MikkelSchubert/adapterremoval
AdapterRemoval --file1 $FILE1 --file2 $FILE2 --threads 2 \
--basename $sra --trimns --trimqualities --collapse \
--qualitybase 33 --gzip --qualitybase-output 33 --minquality 5

#Rename trimmed files
mv ${sra}.pair1.truncated.gz ${sra}_1.truncated.gz
mv ${sra}.pair2.truncated.gz ${sra}_2.truncated.gz

# QC individual files
fastqc ${sra}_1.truncated.gz
fastqc ${sra}_2.truncated.gz

unzip ${sra}_1.truncated_fastqc.zip
unzip ${sra}_2.truncated_fastqc.zip

grep "Per base sequence quality" ${sra}_1.truncated_fastqc/summary.txt >> qc_log.txt
grep "Adapter Content" ${sra}_1.truncated_fastqc/summary.txt >> qc_log.txt

grep "Per base sequence quality" ${sra}_2.truncated_fastqc/summary.txt >> qc_log.txt
grep "Adapter Content" ${sra}_2.truncated_fastqc/summary.txt >> qc_log.txt

rm -r ${sra}_1.truncated_fastqc
rm -r ${sra}_2.truncated_fastqc

else

echo "OK"

fi


# END
