#!/bin/bash
#SBATCH --job-name=S3a_genus-species
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=

####notes and usage####
#
##Notes##
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/)
#Edit the step3a_mapping.sh script:edit/add variables "genus_species", "A", "t", "n", "array", "cpus_per_task", "mail_user", "mail_type"
#
#Example: 
#
#for variables:
#genus_species=Marmota-marmota-marmota
#A=fnrdewoody
#t=10-00:00:00
#n=64
#array=1-7%4 (=In case of 7 samples. In case of 25 samples, the maximum number should be 25 (excluding the first header line). "1-7%4" means 4 tasks will be run in parallel out of the 7 tasks.)
#cpus_per_task=16 (=If you want to assign 16x4 cpus in total. In case you want to use 128 cpus in total, this number should be increased to 32.)
#mail_user=jeon96@purdue.edu
#
#usage:
#sbatch $CLUSTER_SCRATCH/theta/step3a_mapping.sh
#
####end usage and notes####

################################
#Enter undefined variable below#
################################

genus_species=
A=
t=
n=
array=
cpus_per_task=
mail_user=

########################
#DO NOT EDIT BELOW CODE#
########################

# first job to preprocess reference genome and create cleaned_sralist 
job_id1=$(sbatch --job-name=S3_1_$genus_species -A $A -t $t -n $n --mail-user=$mail_user --mail-type=FAIL $CLUSTER_SCRATCH/theta/array_step3_1.sh $genus_species | sed 's/Submitted batch job //')
echo "$job_id1 started, check .err and .out files of it for progress"

# second job to map sra fastq files and create final .bam files 
job_id2=$(sbatch --job-name=S3_2_$genus_species -A $A -t $t --array=$array --cpus-per-task=$cpus_per_task --mail-user=$mail_user --mail-type=FAIL --dependency=afterok:$job_id1 $CLUSTER_SCRATCH/theta/array_step3_2.sh $genus_species $cpus_per_task | sed 's/Submitted batch job //')
echo "$job_id2 started, check .err and .out files of it for progress"

# concatenate .err file and .out files
sbatch --job-name=S3c_$genus_species -A $A -t $t -n 1 --mail-user=$mail_user --mail-type=FAIL,END --dependency=afterany:$job_id1:$job_id2 $CLUSTER_SCRATCH/theta/array_step3_cat.sh $genus_species $job_id1 $job_id2

# END
