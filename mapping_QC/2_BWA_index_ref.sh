#!/bin/sh

#PBS -N BWA_index_ref
#PBS -q standby
#PBS -l nodes=1:ppn=1,naccesspolicy=singleuser
#PBS -l walltime=4:00:00
#PBS -m abe

module load bioinfo
module load bwa
module load samtools
module load picard-tools

cd $PBS_O_WORKDIR

gunzip *.fna.gz

mv *fna ref.fa

# index reference and map reads
bwa index -a bwtsw ref.fa

# create reference that reads can be mapped to. Will produce .fai file
samtools faidx ref.fa

PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

# END
