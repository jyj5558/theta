#!/bin/bash
#SBATCH --job-name=S5a_genus-species
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=

####notes and usage####
#
##Notes##
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/)
#Edit the step5a_theta.sh script:edit/add variables "genus_species", "A", "t", "n", "array", "cpus_per_task", "mail_user", "mail_type"
#
#Example: 
#
#for variables:
#genus_species=Marmota-marmota-marmota
#A=fnrdewoody
#t=10-00:00:00
#n=128
#array=1-7%4 (=In case of 7 samples. In case of 25 samples, the maximum number should be 25 (excluding the first header line). "1-7%4" means 4 tasks will be run in parallel out of the 7 tasks.)
#cpus_per_task=32 (=If you want to assign 32x4 cpus in total. In case you want to use 64 cpus in total, this number should be reduced to 16.)
#mail_user=jeon96@purdue.edu
#accession=GCA_014523465.1
#
#usage:
#sbatch $CLUSTER_SCRATCH/theta/step5a_theta.sh
#
####end usage and notes####

################################
#Enter undefined variable below#
################################

genus_species=r
accession=
A=
t=
n=
array=
cpus_per_task=
mail_user=

########################
#DO NOT EDIT BELOW CODE#
########################

# first job to estimate GL, SFS, theta, ROH, etc. 
job_id1=$(sbatch --job-name=S5_1_$genus_species -A $A -t $t -n $n --mail-user=$mail_user --mail-type=FAIL $CLUSTER_SCRATCH/theta/array_step5_1.sh $genus_species $n $accession | sed 's/Submitted batch job //')
echo "$job_id1 started, check .err and .out files of it for progress"

# second job to estimate individual heterozygosity values 
job_id2=$(sbatch --job-name=S5_2_$genus_species -A $A -t $t --array=$array --cpus-per-task=$cpus_per_task --mail-user=$mail_user --mail-type=FAIL --dependency=afterany:$job_id1 $CLUSTER_SCRATCH/theta/array_step5_2.sh $genus_species $cpus_per_task | sed 's/Submitted batch job //')
echo "$job_id2 started, check .err and .out files of it for progress"

# third job to map sra fastq files and create final .bam files 
job_id3=$(sbatch --job-name=S5_3_$genus_species -A $A -t $t -n $n --mail-user=$mail_user --mail-type=FAIL --dependency=afterok:$job_id2 $CLUSTER_SCRATCH/theta/array_step5_3.sh $genus_species | sed 's/Submitted batch job //')
echo "$job_id3 started, check .err and .out files of it for progress"

# concatenate .err file and .out files
sbatch --job-name=S5c_$genus_species -A $A -t $t -n 1 --mail-user=$mail_user --mail-type=FAIL,END --dependency=afterany:$job_id1:$job_id2:$job_id3 $CLUSTER_SCRATCH/theta/array_step5_cat.sh $genus_species $job_id1 $job_id2 $job_id3

# END
