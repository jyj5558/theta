#!/bin/bash

# script was written to be run with SLURM job scheduler

# define variables
genus_species=$1
job_id1=$2
job_id2=$3

# move to the directory
cd /path/to/theta/${genus_species}

# concatenate .err and .out files
for f in `ls S2_1_${genus_species}_${job_id1}.err S2_2_${genus_species}_${job_id2}_[0-9]*.err | sort`
do echo $f; cat "$f"; echo ; done > S2_${genus_species}_${job_id1}_${job_id2}.err

for f in `ls S2_1_${genus_species}_${job_id1}.out S2_2_${genus_species}_${job_id2}_[0-9]*.out | sort`
do echo $f; cat "$f"; echo ; done > S2_${genus_species}_${job_id1}_${job_id2}.out

rm S2_1_${genus_species}_${job_id1}.err S2_2_${genus_species}_${job_id2}_[0-9]*.err
rm S2_1_${genus_species}_${job_id1}.out S2_2_${genus_species}_${job_id2}_[0-9]*.out

# END
