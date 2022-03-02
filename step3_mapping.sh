#!/bin/bash
#SBATCH --job-name=align
#SBATCH -A fnrquail
#SBATCH -t 4-00:00:00 
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mem=50GB


module load bioinfo
module load bwa
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bamtools

#########################################################################
# this script maps read files to the QC reference assembly					#
# the following files to use in downstream analyses are created:	  	#
# *.sam file for each SRR*												#
#########################################################################
####notes and usage####
#
##Notes##
# add target species "genus-species"
#Example:
#genus_species=Marmota-marmota-marmota
#theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/)
#usage:
#/scratch/bell/$USER/theta/step3_mapping.sh
#
#
####end usage and notes####

#designate target species
genus_species=Tursiops-aduncus

#Make a new directory to house processed alignment files
mkdir /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/

#move to cleaned fastq files in preparation for alignment
cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/

#Capture all cleaned fastq files with variable
for i in `ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq`
do


# index reference 
bwa index -a bwtsw ../../*_ref/ref.fa

#and  index reference sequence in preparation for step4
samtools faidx ../../*_ref/ref.fa

#perform alignment using twenty CPUs and bwa mem algorithm
rm ../aligned/*sam
bwa mem -t 20 -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ../../*_ref/ref.fa  ${i}_1_val_1.fq ${i}_2_val_2.fq > ../aligned/${i}.sam

#Move to the directory containing the alignment files
cd /scratch/bell/dewoody/theta/${genus_species}/sra/aligned/

#Validated *sam files should produce a summary output file that contains no errors
PicardCommandLine ValidateSamFile I=${i}.sam MODE=SUMMARY O=${i}.sam.txt

# check for errors in sam file
error=$(grep "No errors found" ${i}_samfile.txt)
[[ ! -z "$error" ]] && echo "$i samfile not OK" || echo "$i samfile OK"

#Sort sam file based upon read coordinate
PicardCommandLine SortSam INPUT=${i}.sam OUTPUT=${i}_sorted.bam SORT_ORDER=coordinate

#Mark PCR duplicates without removing them
PicardCommandLine MarkDuplicates INPUT=${i}_sorted.bam OUTPUT=./${i}_marked.bam METRICS_FILE=${i}_metrics.txt
PicardCommandLine BuildBamIndex INPUT=./${i}_marked.bam

rm -rf ../../*_ref/ref.dict
PicardCommandLine CreateSequenceDictionary reference=../../*_ref/ref.fa output=../../*_ref/ref.dict

# local realignment of reads
rm -rf forIndelRealigner.intervals
GenomeAnalysisTK -nt 20 -T RealignerTargetCreator -R ../../*_ref/ref.fa -I ${i}_marked.bam -o forIndelRealigner.intervals

#Realign with established intervals
GenomeAnalysisTK -T IndelRealigner -R ../../*ref/ref.fa -I ${i}_marked.bam -targetIntervals forIndelRealigner.intervals -o ../final_bams/${i}.bam

cd /scratch/bell/dewoody/theta/${genus_species}/sra/aligned/final_bams/

#Get some summar stats on bam files
rm *.txt
samtools flagstat ./${i}.bam > ./${i}_mapping.txt
samtools depth -a ./${i}.bam | awk '{c++;s+=$3}END{print s/c}' > ./${i}_stats.txt
samtools depth -a ./${i}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ./${i}_depth.txt

cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/
done

# END
