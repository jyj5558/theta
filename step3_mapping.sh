#!/bin/bash
#SBATCH --job-name=S3_Genus-species
#SBATCH -A fnrquail
#SBATCH -t 10-00:00:00 
#SBATCH -N 1 
#SBATCH -n 20
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mem=20GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="your email address (e.g.jeon96@purdue.edu) without quotation marks"


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

##########################
#designate target species#
##########################

genus_species=

########################
#Do not edit code below#
########################

#Make a new directory to house processed alignment files
mkdir /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/

#Move to reference

#Index/create dict files
cd /scratch/bell/dewoody/theta/${genus_species}/*_ref/

#Create dictionary for realignment
PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

#move to cleaned fastq files in preparation for alignment
cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/

# index reference 
bwa index -a bwtsw ../../*_ref/ref.fa

#and  index reference sequence in preparation for step4
samtools faidx ../../*_ref/ref.fa

#Capture all cleaned fastq files with variable
for i in `ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq`
do

#perform alignment using twenty CPUs and bwa mem algorithm
bwa mem -t 20 -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ../../*_ref/ref.fa  ${i}_1_val_1.fq ${i}_2_val_2.fq > ../aligned/${i}.sam

#Move to the directory containing the alignment files
cd /scratch/bell/dewoody/theta/${genus_species}/sra/aligned/

#Validated *sam files should produce a summary output file that contains no errors
PicardCommandLine ValidateSamFile I=${i}.sam MODE=SUMMARY O=${i}.sam.txt

# check for errors in sam file
noerror=$(grep "No errors found" ${i}.sam.txt)
[[ ! "$noerror" ]] && echo "$i samfile not OK" || echo "$i samfile OK"

#Sort sam file based upon read coordinate
PicardCommandLine SortSam INPUT=${i}.sam OUTPUT=${i}_sorted.bam SORT_ORDER=coordinate

#Mark PCR duplicates without removing them
PicardCommandLine MarkDuplicates INPUT=${i}_sorted.bam OUTPUT=./${i}_marked.bam METRICS_FILE=${i}_metrics.txt
PicardCommandLine BuildBamIndex INPUT=./${i}_marked.bam

# local realignment of reads
GenomeAnalysisTK -nt 20 -T RealignerTargetCreator -R ../../*_ref/ref.fa -I ${i}_marked.bam -o ${i}_forIndelRealigner.intervals

#Realign with established intervals
GenomeAnalysisTK -T IndelRealigner -R ../../*ref/ref.fa -I ${i}_marked.bam -targetIntervals ${i}_forIndelRealigner.intervals -o ../final_bams/${i}.bam

cd /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/

#Get some summar stats on bam files
echo "{i} mapping rate is" > ./${i}_mapping.txt
samtools flagstat ./${i}.bam >> ./${i}_mapping.txt
echo "{i} depth is" > ./${i}_depth.txt
samtools depth -a ./${i}.bam | awk '{c++;s+=$3}END{print s/c}' >> ./${i}_depth.txt
echo "{i} breadth is" > ./${i}_breadth.txt
samtools depth -a ./${i}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./${i}_breadth.txt

cd /scratch/bell/dewoody/theta/${genus_species}/sra/cleaned/
done

# END
