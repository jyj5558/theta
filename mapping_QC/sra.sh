#!/bin/bash
#SBATCH --job-name=sra
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo

# cd $SLURM_SUBMIT_DIR
rm -rf *gz

sed -e 's/^M//g' sra.txt > infile2.txt

#### SRA ####

tail -n +2 infile2.txt | while read -r LINE

do

LINE=${LINE%$'\n'}

echo $LINE

acno=${LINE:3}
lng=${#acno}

first6=${LINE:0:6}

if [ $lng -lt 7 ]
then

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${LINE}/*.fastq.gz

elif [ $lng -eq 7 ]
then

subdir=00${LINE:(-1)}

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${subdir}/${LINE}/*.fastq.gz

elif [ $lng -eq 8 ]
then

subdir=0${LINE:(-2)}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${subdir}/${LINE}/*.fastq.gz

elif [ $lng -eq 9 ]
then

subdir=${LINE:(-3)}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${subdir}/${LINE}/*.fastq.gz

fi

done

# END
