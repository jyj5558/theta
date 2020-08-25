#!/bin/sh

#PBS -N picard
#PBS -q fnrdewoody
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=300:00:00
#PBS -m abe

module load bioinfo
module load samtools
module load picard-tools/2.9.0

cd $PBS_O_WORKDIR

# check coverage
PicardCommandLine CollectWgsMetrics I=marked.bam O=coverage_marked.txt R=ref.fa 


# END
