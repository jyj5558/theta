#!/bin/bash
#SBATCH --job-name=realign_auto
#SBATCH -A fnrdewoody
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bioawk

# cd $SLURM_SUBMIT_DIR

rm -rf *marked*

# this step assures there are no bam errors (hopefully)
#PicardCommandLine SamFormatConverter I=out.bam O=out.sam

#PicardCommandLine SamFormatConverter I=out.sam O=corrected.bam

# check bam file
#PicardCommandLine ValidateSamFile       I=out.bam MODE=SUMMARY

# sort bam file
PicardCommandLine  SortSam I=merged_sorted.bam O=sorted.bam SORT_ORDER=coordinate

# marking PCR duplicated reads without removing them
PicardCommandLine MarkDuplicates INPUT=sorted.bam OUTPUT=marked.bam M=metrics.txt

# check coverage
PicardCommandLine CollectWgsMetrics I=marked.bam O=coverage_marked.txt R=ref.fa 

PicardCommandLine BuildBamIndex INPUT=marked.bam

# create reference that reads can be mapped to. Will produce .fai file
samtools faidx ref.fa

rm -rf ref.dict

PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

#GATK was used for local realignment of reads
GenomeAnalysisTK -nt 1 -T RealignerTargetCreator -R ref.fa -I marked.bam \
-o forIndelRealigner.intervals

GenomeAnalysisTK -T IndelRealigner -R ref.fa -I marked.bam \
-targetIntervals forIndelRealigner.intervals -o realigned_reads.bam

##############################
# make bed and bam file #

# remove files that are not used
rm -rf round*
rm -rf *100k*
rm -rf *.e*
rm -rf *.o*
rm -rf autosomes*
rm -rf tmp

mkdir tmp

# subset ref to required legnth
bioawk -c fastx '{ if(length($seq) > 10000) { print ">"$name; print $seq }}' ref.fa > ref_10k.fa

# index ref
samtools faidx ref_10k.fa

# make scaffold list
cut -f1 ref_10k.fa.fai > scaffold

# make list with sex scafs 
cat mtDNA/satsuma_summary.chained.out chrZ10k/satsuma_summary.chained.out \
chrW10k/satsuma_summary.chained.out > list10k.txt

# take the 4th column
cut -f4 list10k.txt > list10k

# find the unique scaffolds
sort -u list10k > unique_10k.txt

# find differences
sort scaffold unique_10k.txt | uniq -u > autosomes.txt

# split to make sure all autosome scaffolds are represented
rm -rf x*

split -l 5000 autosomes.txt 
# make list of subfiles
ls x* > auto.txt

cat auto.txt | while read -r LINE
do
grep -f <(sed 's/.*/\^&\\>/' $LINE) ref_10k.fa.fai >> autosomes.fa.fai
done

# make bed file 
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' autosomes.fa.fai > autosomes.bed

# make sure that autosomes.txt and autosomes.fa.fai has the same number of lines
autotxt=$(cat autosomes.txt | wc -l)
autofai=$(cat autosomes.fa.fai | wc -l)

#if [ "$autotxt" == "$autofai" ]
#then

# make bed file 
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' autosomes.fa.fai > autosomes.bed

# change RG in bam
PicardCommandLine AddOrReplaceReadGroups \
      I=realigned_reads.bam \
      O=RG.bam \
      RGID=group1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample1
#      "@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1"

# https://broadinstitute.github.io/picard/explain-flags.html
samtools view -bh -f 3 -L autosomes.bed RG.bam > autosomes.bam

# sort bam file
samtools sort -T tmp autosomes.bam -o autosomes_sorted.bam

# convert to sam
samtools view -h -o autosomes.sam autosomes.bam

# check that the new bam is not corrupted
PicardCommandLine ValidateSamFile I=autosomes.sam IGNORE_WARNINGS=true \
MODE=SUMMARY O=autosomes_out.txt TMP_DIR=`pwd`/tmp

rm -rf autosomes_${LINE}.sam
rm -rf autosomes_${LINE}_sorted*

#else
#  echo "autosome files differ"
#fi  

# END
