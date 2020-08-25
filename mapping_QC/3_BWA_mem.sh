#!/bin/bash
#SBATCH --job-name=BWA
#SBATCH -A fnrdewoody
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 5 

module load bioinfo
module load bwa

# cd $SLURM_SUBMIT_DIR

echo $PWD > pwd

sed -i "s/trimgalore/ref.fa ./g" pwd

pwd=$(<pwd)

cp $pwd

bwa index -a bwtsw ref.fa

# fastq
bwa mem -t 5 -M -R "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1" \
ref.fa read1_val_1.fq.gz read2_val_2.fq.gz > mapped.sam

rm -rf pwd

# END
