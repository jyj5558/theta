#!/bin/bash
#SBATCH --job-name=renameFA
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools

# cd $SLURM_SUBMIT_DIR

# rename bam header
samtools view -H test.bam | sed 's/\.1__.*/.1/' | samtools reheader - test.bam > test2.bam

# change mappability file
cut -f1 mappability.bedgraph > scafs
sed -i 's/\.1__.*/.1/g' scafs
paste scafs mappability.bedgraph > mappability.test
cat mappability.test |awk '{print $1,$3,$4,$5}'|sed 's/ /\t/g' > mappability.test2