#!/bin/bash
#SBATCH --job-name=angsdFish
#SBATCH -A fnrquail
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load gsl
module load zlib
module load gcc
module load samtools
module load bioawk

# cd $SLURM_SUBMIT_DIR

export PATH=/home/abruenic/angsd/:$PATH
export PATH=/home/abruenic/angsd/misc/:$PATH

# index bam and ref
samtools index in.bam
samtools faidx ref.fa

# convert bed file to angsd format
awk '{print $1"\t"$2+1"\t"$3}' ok.bed > angsd.file

# index file
angsd sites index angsd.file

# estimate GL
angsd -i in.bam -sites angsd.file -anc ref.fa -dosaf 1 -gl 1 \
-minMapQ 30 -minQ 20 -rf chrs.txt -out out

# obtain ML estimate of SFS using the realSFS
realSFS out.saf.idx > out.sfs

# calculate theta for each site
realSFS saf2theta out.saf.idx -sfs out.sfs -outname out

# estimate Tajimas D
thetaStat do_stat out.thetas.idx

thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz

# column 4 has Wattersons
awk '{print $4}' theta.thetasWindow.gz.pestPG > Watterson

# get mean
meanW=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' Watterson)

# get SD
sdW=$(awk '{delta = $1 - avg; avg += delta / NR; \
mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' Watterson)

# print to file
echo -e "$PWD\t $meanW\t $sdW" \
>> /scratch/snyder/a/abruenic/pupfish/Wattersons_theta.txt

# END