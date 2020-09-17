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
module load BEDTools

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script is identifying repeats, mappability and short scaffolds 
# the following files to use in downstream analyses are created:	  	
# ref.fa (reference file with scaffolds>100kb)							
# ok.bed (regions to analyze in angsd etc)									
# check for repeatmasker file on NCBI to skip that step (*_rm.out.gz )	
# https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/			
#########################################################################

# make sure you have genmap installed and in your path
# https://github.com/cpockrandt/genmap
export PATH=/home/abruenic/genmap-build/bin:$PATH

# download reference from NCBI (Black circulated script for this)
# wget "LINK TO REF" .

# insert ref name here
REF=$"****.fna"

# unzip ref and change name
# keep a copy of the original reference. We to be use which version we are using.
gunzip ${REF}.gz
cp ${REF} original.fa

# sort by lenght
sortbyname.sh in=original.fa out=sorted.fa length descending

# sort ref and rename scaffolds
bioawk -c fastx '{ print ">scaffold-" ++i" "length($seq)"\n"$seq }' \
< sorted.fa > ref.fa

# replace gap with _ in scaffold names
sed -i 's/ /_/g' ref.fa

# make table with old and new scaffold names
paste <(grep ">" sorted.fa) <(grep ">" ref.fa) | sed 's/>//g' \
> scaffold_names.txt

# identify repeats in ref
RepeatMasker -qq -species mammal ref.fa

# build an index of the fasta file(s) whose mappability you want to compute
genmap index -F ref.fa -I index -S 30

# compute mappability
# k = kmer of 100bp
# E = # two mismatches
genmap map -K 100 -E 2 -I index -O mappability -t -w -bg

# make bed file from repeatmasker
cat ref.fa.out|tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' \
> repeats.bed

# sort bed 
sort -V -k 1,3 "repeats.bed" | sortBed | tee test2.bed | sort -c -k1,1 -k2,2n || \
true > repeats_sorted.bed

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

# remove scaffolds shorter than 100kb
bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' \
ref.fa > ref_100k.fa

# index ref
samtools faidx ref_100k.fa

# make list with the >100kb scaffolds
cut -f1 ref_100k.fa.fai |sort|uniq > chrs.txt

# only include chr in merged.bed if they are in chrs.txt
grep -f chrs.txt merged.bed > ok.bed

# remove excess files
rm -rf original.fa
rm -rf sorted.fa
rm -rf merged.bed
rm -rf filter.bed
rm -rf map.bed
rm -rf ref_sorted.genome
rm -rf ref.genome
rm -rf repeats_sorted.bed
rm -rf repeats.bed
rm -rf scaffold_names.txt
rm -rf sorted.fa

# END