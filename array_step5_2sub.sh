#!/bin/bash

#######################################################################
# this script identifies/filteres for usable sites and estimates theta#
# the following result files are created:	#
# Het (=individual heterozygosity; will be created by array job)
#######################################################################

########################
#Do not edit code below#
########################

# load modules
module --force purge
module load biocontainers
module load angsd
module load r

# define variables
genus_species=$1
bam=$2
cpus_per_task=$3

echo "The task $SLURM_ARRAY_TASK_ID handles $bam"

# path to parent directory of genus-species
PD=/scratch/bell/dewoody/theta/${genus_species}/

# just some absolute  paths to shorten commands
FINAL=/scratch/bell/dewoody/theta/${genus_species}/sra/final_bams
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta

# move to directory which will house files
cd ${THETA}

# calculate individual heterozygosity of this ${bam}

OUTDIR='HET'

echo "${bam} heterozygosity estimation started"

angsd -i ${FINAL}/${bam}.bam -ref ${PD}/*_ref/ref.fa -anc ${PD}/*_ref/ref.fa -dosaf 1 -rf ${PD}/*_ref/chrs.txt -sites ./angsd.file \
-minMapQ 30 -minQ 30 -P ${cpus_per_task} -out ${OUTDIR}/${bam} -only_proper_pairs 1 -baq 2 \
-GL 2 -doMajorMinor 1 -doCounts 1 -setMinDepthInd 5 -uniqueOnly 1 -remove_bads 1 

realSFS -P ${cpus_per_task} -fold 1 ${OUTDIR}/${bam}.saf.idx > ${OUTDIR}/${bam}_est.ml

cd ${OUTDIR}
Rscript -e 'args<-commandArgs(TRUE); bam<-args[1]; a<-scan(paste(bam,"est.ml", sep="_")); a[2]/sum(a)' ${bam} >>  ../Het

echo "${bam} heterozygosity estimation done"
cd ../

# END
