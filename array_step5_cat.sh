#!/bin/bash

# define variables
genus_species=$1
job_id1=$2
job_id2=$3
job_id3=$4

# move to the directory
cd /scratch/bell/dewoody/theta/${genus_species}

# concatenate .err and .out files
for f in `ls S5_1_${genus_species}_${job_id1}.err S5_2_${genus_species}_${job_id2}_[0-9]*.err S5_3_${genus_species}_${job_id3}.err | sort`
do echo $f; cat "$f"; echo ; done > S5_${genus_species}_${job_id1}_${job_id2}_${job_id3}.err

for f in `ls S5_1_${genus_species}_${job_id1}.out S5_2_${genus_species}_${job_id2}_[0-9]*.out S5_3_${genus_species}_${job_id3}.out | sort`
do echo $f; cat "$f"; echo ; done > S5_${genus_species}_${job_id1}_${job_id2}_${job_id3}.out

rm S5_1_${genus_species}_${job_id1}.err S5_2_${genus_species}_${job_id2}_[0-9]*.err S5_3_${genus_species}_${job_id3}.err
rm S5_1_${genus_species}_${job_id1}.out S5_2_${genus_species}_${job_id2}_[0-9]*.out S5_3_${genus_species}_${job_id3}.out

# END
