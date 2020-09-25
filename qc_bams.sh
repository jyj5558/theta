#!/bin/bash
#SBATCH --job-name=qcBam
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

samtools index ${LINE}

done

# global depth (read count across all samples)
angsd -bam bam.filelist -doDepth 1 -out strict -doCounts 1 -minMapQ 30 \
-minQ 20 -maxDepth 600

# find threshold excluding the sites with 1% lowest and 1% highest global depth 
Rscript /scratch/snyder/a/abruenic/scripts/percentiles.R

q1=$(head -n 1 strict.percentile)
q2=$(tail -n 1 strict.percentile)


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