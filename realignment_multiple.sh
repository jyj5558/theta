#!/bin/bash
#SBATCH --job-name=realignment
#SBATCH -A fnrquail
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bamtools

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script check for errors, mark duplicates and realign bam files
#########################################################################

# make sure you have the SRAs for each sample in the sample DIR

cat samples.txt | while read -r SAMPLE

do

cd $SAMPLE

# make tmp directory for the files (remove potential "old" tmps)
rm -rf tmp
mkdir tmp

# make list of the samfiles
ls *.sam > samfiles.txt

sed -i "s/.sam//g" samfiles.txt

cat samfiles.txt | while read -r LINE

do

# validate the SAM file should produce a validate_output.txt file 
# that says there are no errors.
PicardCommandLine ValidateSamFile I=${LINE}.sam MODE=SUMMARY O=${LINE}_samfile.txt
PicardCommandLine SortSam SORT_ORDER=coordinate INPUT=${LINE}.sam OUTPUT=${LINE}.bam \
TMP_DIR=`pwd`/tmp VALIDATION_STRINGENCY=LENIENT

# check for errors in sam file
error=$(grep "No errors found" ${LINE}_samfile.txt)

[[ ! -z "$error" ]] && echo "$LINE samfile not OK" || echo "$LINE samfile OK"

# perform realignment etc on each individual
# marking PCR duplicated reads without removing them
PicardCommandLine MarkDuplicates INPUT=${LINE}.bam OUTPUT=${LINE}_marked.bam \
M=metrics.txt

PicardCommandLine BuildBamIndex INPUT=${LINE}_marked.bam

# index reference
samtools faidx ref.fa

rm -rf ref.dict
PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

rm -rf forIndelRealigner.intervals
# local realignment of reads
GenomeAnalysisTK -nt 1 -T RealignerTargetCreator -R ref.fa -I ${LINE}_marked.bam \
-o forIndelRealigner.intervals

GenomeAnalysisTK -T IndelRealigner -R ref.fa -I ${LINE}_marked.bam \
-targetIntervals forIndelRealigner.intervals -o ${LINE}_realigned_reads.bam

done

cd ..

done

# remove excess files
rm -rf *marked.ba*
rm -rf out.ba*
rm -rf *.sam
rm -rf forIndelRealigner.intervals

# END