#!/bin/bash
#SBATCH --job-name=S2a_genus-species	
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=

####notes and usage####
#
##Notes##
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/)
#Edit the step3a_mapping.sh script:edit/add variables "genus_species", "A", "t", "n", "array", "cpus_per_task", "mail_user", "mail_type"
#script was written to be run with SLURM job scheduler
#
#Example: 
#
#for variables:
#genus_species=Marmota-marmota-marmota
#A=
#t=10-00:00:00
#n=128
#array=1-7%4 (=In case of 7 samples. In case of 25 samples, the maximum number should be 25 (excluding the first header line). "1-7%4" means 4 tasks will be run in parallel out of the 7 tasks.)
#cpus_per_task=32 (=If you want to assign 32x4 cpus in total. In case you want to use 64 cpus in total, this number should be reduced to 16.)
#mail_user=
#mail_type=BEGIN,FAIL,END
#
#usage:
#sbatch $CLUSTER_SCRATCH/theta/step2a_SRA_download.sh
#
####end usage and notes####

################################
#Enter undefined variable below#
################################

genus_species=Ovis-canadensis
A=
t=
n=
array=
cpus_per_task=
mail_user=

########################
#DO NOT EDIT BELOW CODE#
########################

# first job to create sample lists and directories 
job_id1=$(sbatch --job-name=S2_1_$genus_species -A $A -t $t -n $n --mail-user=$mail_user --mail-type=FAIL $CLUSTER_SCRATCH/theta/array_step2_1.sh $genus_species | sed 's/Submitted batch job //')
echo "$job_id1 submitted, check .err and .out files of it for progress"

# second job to download and clean sra files
job_id2=$(sbatch --job-name=S2_2_$genus_species -A $A -t $t --array=$array --cpus-per-task=$cpus_per_task --mail-user=$mail_user --mail-type=FAIL --dependency=afterany:$job_id1 $CLUSTER_SCRATCH/theta/array_step2_2.sh $genus_species $cpus_per_task | sed 's/Submitted batch job //')
echo "$job_id2 submitted, check .err and .out files of it for progress"

# concatenate .err file and .out files
sbatch --job-name=S2c_$genus_species -A $A -t $t -n 1 --mail-user=$mail_user --mail-type=FAIL,END --dependency=afterany:$job_id1:$job_id2 $CLUSTER_SCRATCH/theta/array_step2_cat.sh $genus_species $job_id1 $job_id2

# END
