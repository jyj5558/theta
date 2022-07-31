#!/bin/bash 
#SBATCH -e %x_%A_%a.err
#SBATCH -o %x_%A_%a.out

############################################################################
# this script will create:
# raw and cleaned .fastq files for each SRA using array_step2_2sub.sh script
############################################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# defnie argument variables
genus_species=$1
cpus_per_task=$2

# capture all sra accession numbers with variable
sra=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/bell/dewoody/theta/${genus_species}/sra/raw/sralist)

# download each sra files by running array job
srun $CLUSTER_SCRATCH/theta/array_step2_2sub.sh ${genus_species} ${sra} ${cpus_per_task}
