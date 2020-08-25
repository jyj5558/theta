#!/bin/bash
#SBATCH --job-name=realign_auto
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load samtools
module load picard-tools/2.9.0
module load GATK/3.8.1
module load bioawk
module load bwa

# cd $SLURM_SUBMIT_DIR
rm -rf auto*
rm -rf tmp
mkdir tmp

bwa index -a bwtsw ref.fa

samtools faidx ref.fa

ls *.sam > samfiles.txt

sed -i "s/.sam//g" samfiles.txt

cat samfiles.txt | while read -r LINE

do

# validate the SAM file should produce a validate_output.txt file that says there are no errors.
# samtools view -h -o aln.sam out.bam 
# one section for each SRA is needed. 
# Easiest to copy/paste and use find/replace to insert SRA numbers.
PicardCommandLine ValidateSamFile I=${LINE}.sam MODE=SUMMARY O=${LINE}_samfile.txt
PicardCommandLine SortSam SORT_ORDER=coordinate INPUT=${LINE}.sam OUTPUT=${LINE}.bam \
TMP_DIR=`pwd`/tmp VALIDATION_STRINGENCY=LENIENT

done

# sort bamfiles
ls *.bam > files.bamlist

sed "s/.bam//g" files.bamlist > bamfiles

cat bamfiles| while read -r LINE

do

samtools sort -o ${LINE}_sorted.bam ${LINE}.bam

done

# make list of sorted bams
ls *_sorted.bam > sorted.bamlist

samtools merge -b sorted.bamlist merged.bam

samtools sort -o merged_sorted.bam merged.bam

samtools flagstat merged_sorted.bam

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


# END
