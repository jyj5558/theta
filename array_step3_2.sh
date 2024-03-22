#!/bin/bash
#SBATCH -e %x_%A_%a.err
#SBATCH -o %x_%A_%a.out

##########################################################
# this script will create:
# .bam files for each SRA using array_step3_2sub.sh script
# script was written to be run with SLURM job scheduler
##########################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# define variables
genus_species=$1
cpus_per_task=$2
sra=$(sed -n "$SLURM_ARRAY_TASK_ID"p /path/to/theta/${genus_species}/sra/cleaned/cleaned_sralist)

# mapping sra reads to the reference genome by running array job
srun $CLUSTER_SCRATCH/theta/array_step3_2sub.sh ${genus_species} ${sra} ${cpus_per_task}

# END
