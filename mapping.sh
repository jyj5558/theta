#!/bin/bash
#SBATCH --job-name=BWA
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 5 

module load bioinfo
module load bwa

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script maps read files to the reference assembly					#
# the following files to use in downstream analyses are created:	  	#
# *.sam file for each SRR*												#
#########################################################################

# index reference and map reads
bwa index -a bwtsw ref.fa

cat fastq.txt | while read -r LINE

do

if [ -d /${LINE} ]
then
# fastq
bwa mem -t 5 -M -R "@RG\tID:group1\tSM:${LINE}\tPL:illumina\tLB:lib1\tPU:unit1" \
ref.fa ${LINE}/${LINE}_1_val_1.fq.gz ${LINE}/${LINE}_2_val_2.fq.gz > ${LINE}.sam

else 
bwa mem -t 5 -M -R "@RG\tID:group1\tSM:${LINE}\tPL:illumina\tLB:lib1\tPU:unit1" \
ref.fa ${LINE}_1.fastq.gz ${LINE}_2.fastq.gz > ${LINE}.sam

fi
done

# END
