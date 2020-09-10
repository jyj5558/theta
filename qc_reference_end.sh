#!/bin/bash
#SBATCH --job-name=qcRef
#SBATCH -A fnrfish
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load bioawk
module load seqtk
module load RepeatMasker/4.0.7
module load samtools
module load cmake/3.9.4
module load BBMap

# cd $SLURM_SUBMIT_DIR

# make sure you have genmap installed and in your path
# https://github.com/cpockrandt/genmap
export PATH=/home/abruenic/genmap-build/bin:$PATH

# download reference
# wget "LINK TO REF" .


# identify repeats in ref
RepeatMasker ref.fa

# build an index of the fasta file(s) whose mappability you want to compute
genmap index -F ref.fa -I index -S 30

# compute mappability
# k = kmer of 100bp
# E = # two mismatches
genmap map -K 100 -E 2 -I index -O mappability -t -w -bg

# make bed file from repeatmasker
cat ref.fa.out|tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' \
> repeats.bed

# index ref
samtools faidx ref.fa

# make ref.genome
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ref.fa.fai > ref.genome

# sort genome file
sort -k1 ref.genome > ref_sorted.genome

# find nonrepeat regions
bedtools complement -i repeats_sorted.bed -g ref_sorted.genome > nonrepeat.bed

# clean mappability file, remove sites with <1 mappability
awk '$4 == 1' mappability.bedgraph > map.bed

awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' map.bed > mappability.bed

# combine bed files for mappability and repeatmasker
cat mappability.bed nonrepeat.bed > filter.bed

# sort combined file -- by chr then site
bedtools sort -i filter.bed > filter_sorted.bed

# merge overlapping regions
bedtools merge -i filter_sorted.bed > merged.bed

# make list with the >100kb scaffolds
cut -f1 ref.fa.fai |sort|uniq >chrs.txt

# only include chr in merged.bed if they are in chrs.txt
grep -f chrs.txt merged.bed > ok.bed

# END