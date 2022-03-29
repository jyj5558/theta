#!/bin/bash
#SBATCH --job-name=qcRef
#SBATCH -A fnrfdewoody
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mem=50G

module load bioinfo
module load bioawk
module load seqtk
module load RepeatMasker/4.0.7
module load samtools
module load cmake/3.9.4
module load BEDTools
module load BBMap
module load r

# cd $SLURM_SUBMIT_DIR

#########################################################################
# this script is identifying repeats, mappability and short scaffolds 
# the following files to use in downstream analyses are created:	  	
# ref.fa (reference file with scaffolds>100kb)							
# ok.bed (regions to analyze in angsd etc)									
# check for repeatmasker file on NCBI to skip that step (*_rm.out.gz )
# https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#########################################################################

# make sure you have genmap installed and in your path
# https://github.com/cpockrandt/genmap
export PATH=/home/abruenic/genmap-build/bin:$PATH

# download reference from NCBI
# see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
wget "LINK TO REF" .
gunzip *fna.gz

# keep a copy of the original reference.
cp *.fna original.fa

# download repeatmasker from NCBI
# see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
wget "LINK TO REPEATMASKER" .
gunzip *_rm.out.gz

# keep a copy of the original repeatmasker. 
cp *_rm.out rm.out

# sort ref by lenght
sortbyname.sh in=original.fa out=sorted.fa length descending

# sort ref and rename scaffolds
bioawk -c fastx '{ print ">scaffold-" ++i" "length($seq)"\n"$seq }' \
< sorted.fa > ref.fa

# replace gap and - with _ in scaffold names
sed -i 's/ /_/g' ref.fa
sed -i 's/-/_/g' ref.fa

samtools faidx ref.fa

# make table with old and new scaffold names
paste <(grep ">" sorted.fa) <(grep ">" ref.fa) | sed 's/>//g' \
> scaffold_names.txt

# make file with ID and scaffold identifier
awk '{ print $1, $NF }' scaffold_names.txt > ID.txt

FILE1=$"rm.out"

if [ -s $FILE1 ]
then
	# change scaffold names to fit ref.fa for rm.out 
	tail -n +4 rm.out > rm.body
	head -n 3 rm.out > rm.header
	sed -i 's/\*//g' rm.body 
	# number of columns
	Rscript /scratch/snyder/a/abruenic/scripts/repeatmasker_names.R
	sed -i 's/"//g' rm_edited.body 
	cat rm.header rm_edited.body > rm.out
	# make bed file from repeatmasker
	cat rm.out |tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' \
	> repeats.bed
else	
	# if no repeatmasker file is available run RepeatMasker
	RepeatMasker -qq -species mammal ref.fa 
	# make bed file from NCBI repeatmasker
	cat repeatmasker.fa.out|tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' \
	> repeats.bed
fi

# build an index of the fasta file(s) whose mappability you want to compute
rm -rf index
genmap index -F ref.fa -I index -S 50

# compute mappability
# k = kmer of 100bp
# E = # two mismatches
genmap map -K 100 -E 2 -I index -O mappability -t -w -bg

# sort bed 
sortBed -i repeats.bed > repeats_sorted.bed

# make ref.genome
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ref.fa.fai > ref.genome

# sort genome file
awk '{print $1, $2, $2}' ref.genome > ref2.genome
sed -i 's/ /\t/g' ref2.genome
sortBed -i ref2.genome > ref3.genome
awk '{print $1, $2 }' ref3.genome > ref_sorted.genome
sed -i 's/ /\t/g' ref_sorted.genome

# find nonrepeat regions
bedtools complement -i repeats_sorted.bed -g ref_sorted.genome > nonrepeat.bed

# clean mappability file, remove sites with <1 mappability
awk '$4 == 1' mappability.bedgraph > map.bed

awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' map.bed > mappability.bed

# sort mappability 
sortBed -i mappability.bed > mappability2.bed
sed -i 's/ /\t/g' mappability2.bed

# only include sites that are nonrepeats and have mappability ==1
bedtools subtract -a mappability2.bed -b repeats_sorted.bed > map_nonreapeat.bed

# sort file -- by chr then site
bedtools sort -i map_nonreapeat.bed > filter_sorted.bed

# merge overlapping regions
bedtools merge -i filter_sorted.bed > merged.bed

# remove scaffolds shorter than 100kb
bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' \
ref.fa > ref_100k.fa

# index
samtools faidx ref_100k.fa

# make list with the >10kb scaffolds
awk '{ print $1, $2, $2 }' ref_100k.fa.fai > chrs.info

# replace column 2 with zeros
awk '$2="0"' chrs.info > chrs.bed

# make tab delimited
sed -i 's/ /\t/g' chrs.bed

# make chrs.txt
cut -f 1 chrs.bed > chrs.txt

# only include chr in merged.bed if they are in chrs.txt
bedtools intersect -a chrs.bed -b merged.bed > ok.bed	
	
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
rm -rf ref2.genome
rm -rf ref3.genome
# END