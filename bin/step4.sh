#!/bin/bash

module load biocontainers
module load cmake
module load bedtools
module load bbmap
module load r
module load bedops
module load genmap
export PATH=$PATH:~/genmap-build/bin

####notes and usage####
#
# This script downloads reference, repeat, and annotation data and then identifyies repeats, 
# estimates mappability and finds all of the short scaffolds. 
# The output files include: 								
# ok.bed (regions to analyze in angsd etc)		
# map_repeat_summary.txt (summary of ref quality)							
# script was written to be run with SLURM job scheduler
#
####end of notes and usage####

########################
#DO NOT EDIT BELOW CODE#
########################

genus_species=$1
accession=$2
pathway=$3
assembly=$4

#move to species/accession directory

cd $CLUSTER_SCRATCH/theta/${genus_species}/

####assess mappability of reference####
genmap index -F ${accession}_ref/ref_100kb.fa -I index -S 50 # build an index 

# compute mappability, k = kmer of 100bp, E = # two mismatches
mkdir -p mappability
genmap map -K 100 -E 2 -T 64 -I index -O mappability -t -w -bg                

# sort bed 
sortBed -i ${accession}_rm/repeats.bed > ${accession}_rm/repeats_sorted.bed 

# make ref.genome
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ${accession}_ref/ref_100kb.fa.fai > ${accession}_ref/ref.genome 

# sort genome file
awk '{print $1, $2, $2}' ${accession}_ref/ref.genome > ${accession}_ref/ref2.genome
sed -i 's/ /\t/g' ${accession}_ref/ref2.genome
sortBed -i ${accession}_ref/ref2.genome > ${accession}_ref/ref3.genome
awk '{print $1, $2 }' ${accession}_ref/ref3.genome > ${accession}_ref/ref_sorted.genome
sed -i 's/ /\t/g' ${accession}_ref/ref_sorted.genome
rm ${accession}_ref/ref.genome
rm ${accession}_ref/ref2.genome
rm ${accession}_ref/ref3.genome

# find nonrepeat regions
bedtools complement -i ${accession}_rm/repeats_sorted.bed -g ${accession}_ref/ref_sorted.genome > ${accession}_rm/nonrepeat.bed

# clean mappability file, remove sites with <1 mappability                                                    
awk '$4 == 1' mappability/ref_100kb.genmap.bedgraph > mappability/map.bed                                           
awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' mappability/map.bed > mappability/mappability.bed

# sort mappability 
sortBed -i mappability/mappability.bed > mappability/mappability2.bed
sed -i 's/ /\t/g' mappability/mappability2.bed

# only include sites that are nonrepeats and have mappability ==1
bedtools subtract -a mappability/mappability2.bed -b ${accession}_rm/repeats_sorted.bed > mappability/map_nonreapeat.bed

# sort file -- by chr then site
bedtools sort -i mappability/map_nonreapeat.bed > mappability/filter_sorted.bed

# merge overlapping regions
bedtools merge -i mappability/filter_sorted.bed > mappability/merged.bed

# make bed file with the 100k and merged.bed (no repeats, mappability =1 sites) from below
awk '{ print $1, $2, $2 }' ./${accession}_ref/ref_100kb.fa.fai > ./${accession}_ref/ref_100kb.info

# replace column 2 with zeros
awk '$2="0"' ./${accession}_ref/ref_100kb.info > ./${accession}_ref/ref_100kb.bed

# make tab delimited
sed -i 's/ /\t/g' ./${accession}_ref/ref_100kb.bed

# only include scaffolds in merged.bed if they are in ref_100kb.bed
bedtools intersect -a ./${accession}_ref/ref_100kb.bed -b ./mappability/merged.bed > ok.bed	

# make chrs.txt
cut -f 1 ${accession}_ref/ref_100kb.bed > ${accession}_ref/chrs.txt
	
# remove excess files
rm mappability/map.bed
rm mappability/filter_sorted.bed
rm mappability/mappability2.bed
rm ${accession}_ref/ref_sorted.genome

#output some QC stats
cd $CLUSTER_SCRATCH/theta/bin/
Rscript qc_reference_stats.R --args $CLUSTER_SCRATCH/theta/${genus_species}/ ${genus_species} ${accession} $CLUSTER_SCRATCH/theta/
cd $CLUSTER_SCRATCH/theta/${genus_species}/

map=$(sed -n '1p' okmap.txt)
norepeats=$(sed -n '1p' norepeat.txt)
okbed=$(sed -n '1p' okbed.txt)

echo -e "${genus_species}\t $accession\t $map\t $norepeats\t $okbed" >> map_repeat_summary.txt

# make list of the bamfiles and index each file in theta directory for step5 in advance
# if you remove some individuals depending on their mapping rate, depth, or breadth, edit this file accordingly before step5
mkdir -p $CLUSTER_SCRATCH/theta/${genus_species}/theta/
cd $CLUSTER_SCRATCH/theta/${genus_species}/theta/
ls $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/*.bam > ./bam.filelist

# END
