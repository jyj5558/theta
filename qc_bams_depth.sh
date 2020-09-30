#!/bin/bash
#SBATCH --job-name=qcBam
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load gsl
module load zlib
module load gcc
module load samtools
module load bioawk
module load r

# cd $SLURM_SUBMIT_DIR

export PATH=/home/abruenic/angsd/:$PATH
export PATH=/home/abruenic/angsd/misc/:$PATH

# index bam and ref
samtools faidx ref.fa

# make list of the bamfiles and index each file
ls *realigned_reads.bam > bam.filelist

cat bam.filelist | while read -r LINE

do

FILE1=$"${LINE}.bai"
if [ -f $FILE1 ]; then
 	echo "OK"
else
	samtools index ${LINE}
fi

done

# global depth (read count across all samples)
angsd -bam bam.filelist -doDepth 1 -out strict -doCounts 1 -minMapQ 30 \
-minQ 20 -dumpCounts 1 -maxdepth 1000 

# find threshold excluding the sites with 1% lowest and 1% highest global depth 
#Rscript /scratch/snyder/a/abruenic/scripts/percentiles.R

