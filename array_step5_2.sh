#!/bin/bash
#SBATCH -e %x_%A_%a.err
#SBATCH -o %x_%A_%a.out

###############################################################################
# this script identifies / filteres for usable sites and estimates theta, etc.#
# the following result files are created: 
# Het for each individual using array_step5_2sub.sh script
###############################################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# define variables
genus_species=$1
cpus_per_task=$2

# first move to the species Theta directory
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta
cd $THETA

# capture all final bam files with variable
ls -1 ../sra/final_bams/*bam | sed 's/\//\'$'\t/g' | cut -f 4| sed 's/.bam//g' > ../sra/final_bams/final_bamlist
bam=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/final_bamlist)

# calculate individual heterozygosity by running array job
srun $CLUSTER_SCRATCH/theta/array_step5_2sub.sh ${genus_species} ${bam} ${cpus_per_task}

# END
