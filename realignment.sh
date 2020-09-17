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

# make tmp directory for the files
mkdir tmp

# make list of the samfiles
ls *.sam > samfiles.txt

sed -i "s/.sam//g" samfiles.txt

cat samfiles.txt | while read -r LINE

do

# validate the SAM file should produce a validate_output.txt file 
# that says there are no errors.
# samtools view -h -o aln.sam out.bam 
# one section for each SRA is needed. 
# Easiest to copy/paste and use find/replace to insert SRA numbers.
PicardCommandLine ValidateSamFile I=${LINE}.sam MODE=SUMMARY O=${LINE}_samfile.txt
PicardCommandLine SortSam SORT_ORDER=coordinate INPUT=${LINE}.sam OUTPUT=${LINE}.bam \
TMP_DIR=`pwd`/tmp VALIDATION_STRINGENCY=LENIENT

done

# check that there are no errors in sam file
cat samfiles.txt | while read -r LINE

do

error=$(grep "No errors found" ${LINE}_samfile.txt)

[[ ! -z "$error" ]] && echo "$LINE samfile not OK" || echo "$LINE samfile OK"

done

# merge bams
# list your SRA numbers in order to merge all bam files into one bam file
ls *.sam > samfiles.txt
sed "s/sam/bam/g" samfiles.txt > files.bamlist 

bamtools merge -list files.bamlist -out out.bam

PicardCommandLine SamFormatConverter I=out.bam O=out.sam

PicardCommandLine SamFormatConverter I=out.sam O=out.bam

# marking PCR duplicated reads without removing them
PicardCommandLine MarkDuplicates INPUT=out.bam OUTPUT=marked.bam M=metrics.txt

# check coverage
PicardCommandLine CollectWgsMetrics I=marked.bam O=coverage_marked.txt R=ref.fa 

PicardCommandLine BuildBamIndex INPUT=marked.bam

# create reference that reads can be mapped to. Will produce .fai file
samtools faidx ref.fa

PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

# local realignment of reads
GenomeAnalysisTK -nt 1 -T RealignerTargetCreator -R ref.fa -I marked.bam \
-o forIndelRealigner.intervals

GenomeAnalysisTK -T IndelRealigner -R ref.fa -I marked.bam \
-targetIntervals forIndelRealigner.intervals -o realigned_reads.bam

# remove excess files
rm -rf marked.bam
rm -rf out.bam
rm -rf *.sam
rm -rf forIndelRealigner.intervals

# END