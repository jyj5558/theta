#!/bin/bash

module purge    
module load biocontainers
module load bwa
module load samtools
module load picard
module load gatk
module load bamtools

####notes and usage####
#
# This script maps read files to the QC reference assembly					#
# the following files to use in downstream analyses are created:	  	#
# *.sam file for each SRR*												#
# script was written to be run with SLURM job scheduler
#
####end of notes and usage####

########################
#DO NOT EDIT BELOW CODE#
########################

genus_species=$1
accession=$2
pathway=$3
assembly=$4
    
#Make a new directory to house processed alignment files
mkdir -p $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/

#Move to reference

#Index/create dict files
cd $CLUSTER_SCRATCH/theta/${genus_species}/${accession}_ref/

#Create dictionary for realignment
PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

#Move to cleaned fastq files in preparation for alignment
cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/cleaned/

#Index reference 
bwa index -a bwtsw ../../${accession}_ref/ref.fa

#Capture all cleaned fastq files with variable
ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq > cleaned_sralist
cat cleaned_sralist
for i in `cat cleaned_sralist`; do

#Perform alignment using whole CPUs and bwa mem algorithm
  bwa mem -t 64 -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ../../${accession}_ref/ref.fa  ${i}_1_val_1.fq ${i}_2_val_2.fq > ../aligned/${i}.sam

#Move to the directory containing the alignment files
  cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/aligned/

#Validated *sam files should produce a summary output file that contains no errors
  PicardCommandLine ValidateSamFile I=${i}.sam MODE=SUMMARY O=${i}.sam.txt

#Check for errors in sam file
  noerror=$(grep "No errors found" ${i}.sam.txt)
  [[ ! "$noerror" ]] && echo "$i samfile not OK" || echo "$i samfile OK"

#Sort sam file based upon read coordinate
  mkdir -p ./tmp/
  PicardCommandLine SortSam TMP_DIR=$CLUSTER_SCRATCH/theta/${genus_species}/sra/aligned/tmp INPUT=${i}.sam OUTPUT=${i}_sorted.bam SORT_ORDER=coordinate

#Mark PCR duplicates without removing them
  PicardCommandLine MarkDuplicates INPUT=${i}_sorted.bam OUTPUT=./${i}_marked.bam METRICS_FILE=${i}_metrics.txt
  PicardCommandLine BuildBamIndex INPUT=./${i}_marked.bam

#Local realignment of reads
  gatk3 -T RealignerTargetCreator -nt 64 -R ../../${accession}_ref/ref.fa -I ${i}_marked.bam -o ${i}_forIndelRealigner.intervals

#Realign with established intervals
  gatk3 -T IndelRealigner -R ../../${accession}_ref/ref.fa -I ${i}_marked.bam -targetIntervals ${i}_forIndelRealigner.intervals -o ../final_bams/${i}.bam
  
  cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/

#Get some summar stats on bam files
  echo "$i mapping rate is" > ./${i}_mapping.txt
  samtools flagstat ./${i}.bam >> ./${i}_mapping.txt
  echo "$i depth is" > ./${i}_depth.txt
  samtools depth -a ./${i}.bam | awk '{c++;s+=$3}END{print s/c}' >> ./${i}_depth.txt
  echo "$i breadth is" > ./${i}_breadth.txt
  samtools depth -a ./${i}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./${i}_breadth.txt

  cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/cleaned/
done

# END