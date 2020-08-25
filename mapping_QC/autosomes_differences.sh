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

rm -rf tmp

mkdir tmp

# subset ref to required legnth
#bioawk -c fastx '{ if(length($seq) > 10000) { print ">"$name; print $seq }}' ref.fa > ref_10k.fa

# index ref
#samtools faidx ref_10k.fa

# make scaffold list
cut -f1 ref_10k.fa.fai > scaffold

# make list with sex scafs 
# this is for the satsuma with no length filter (l = 0)
cat mtDNA/satsuma_summary.chained.out chrZ10k_no_filter/satsuma_summary.chained.out \
chrW10k_no_filter/satsuma_summary.chained.out > newlist10k.txt

# take the 4th column
cut -f4 newlist10k.txt > newlist10k

# find the unique scaffolds
sort -u newlist10k > newunique_10k.txt

# find differences
sort scaffold newunique_10k.txt | uniq -u > autosomes.txt

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

if [ "$autotxt" == "$autofai" ]
then

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

else
  echo $PWD
  echo "autosome files differ"
fi  

# END
