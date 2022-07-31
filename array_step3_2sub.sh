#!/bin/bash

#########################################################################
# this script maps read files to the QC reference assembly					#
# the following files to use in downstream analyses are created:	  	#
# *.bam file for each SRA*												#
#########################################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# load modules
module --force purge
module load bioinfo
module load bwa
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bamtools

# define variables
genus_species=$1
sra=$2
cpus_per_task=$3

echo "The task $SLURM_ARRAY_TASK_ID handles $sra"

# move to the directory containing cleaned_sra files
cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/

# perform alignment using bwa mem algorithm
bwa mem -t ${cpus_per_task} -M -R "@RG\tID:group1\tSM:${sra}\tPL:illumina\tLB:lib1\tPU:unit1" ../../*_ref/ref.fa  ${sra}_1_val_1.fq ${sra}_2_val_2.fq > ../aligned/${sra}.sam

# move to the directory containing the alignment files
cd /scratch/bell/dewoody/theta/${genus_species}/sra/aligned/

# validated *sam files should produce a summary output file that contains no errors
PicardCommandLine ValidateSamFile I=${sra}.sam MODE=SUMMARY O=${sra}.sam.txt

# check for errors in sam file
noerror=$(grep "No errors found" ${sra}.sam.txt)
[[ ! "$noerror" ]] && echo "${sra} samfile not OK" || echo "${sra} samfile OK"

# sort sam file based upon read coordinate
PicardCommandLine SortSam INPUT=${sra}.sam OUTPUT=${sra}_sorted.bam SORT_ORDER=coordinate

# mark PCR duplicates without removing them
PicardCommandLine MarkDuplicates INPUT=${sra}_sorted.bam OUTPUT=./${sra}_marked.bam METRICS_FILE=${sra}_metrics.txt
PicardCommandLine BuildBamIndex INPUT=./${sra}_marked.bam

# local realignment of reads
GenomeAnalysisTK -nt ${cpus_per_task} -T RealignerTargetCreator -R ../../*_ref/ref.fa -I ${sra}_marked.bam -o ${sra}_forIndelRealigner.intervals

# realign with established intervals
GenomeAnalysisTK -T IndelRealigner -R ../../*ref/ref.fa -I ${sra}_marked.bam -targetIntervals ${sra}_forIndelRealigner.intervals -o ../final_bams/${sra}.bam

cd /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/

# get some summar stats on bam files
echo "${sra} mapping rate is" > ./${sra}_mapping.txt
samtools flagstat ./${sra}.bam >> ./${sra}_mapping.txt
echo "${sra} depth is" > ./${sra}_depth.txt
samtools depth -a ./${sra}.bam | awk '{c++;s+=$3}END{print s/c}' >> ./${sra}_depth.txt
echo "${sra} breadth is" > ./${sra}_breadth.txt
samtools depth -a ./${sra}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./${sra}_breadth.txt

# END
