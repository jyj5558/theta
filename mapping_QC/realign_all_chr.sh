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

# cd $SLURM_SUBMIT_DIR

rm -rf auto*

# sort bam file
PicardCommandLine SortSam I=merged_sorted.bam O=sorted.bam SORT_ORDER=coordinate

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
