#!/bin/bash
#SBATCH --job-name=S4_genus-species
#SBATCH -A standby
#SBATCH -t 4:00:00 
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH --mem=100G
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module load bioinfo
module load bioawk
module load cmake/3.9.4
module load BEDTools
module load BBMap
module load R
module load bedops
export PATH=$PATH:~/genmap-build/bin

####notes####
#
# This script downloads reference, repeat, and annotation data and then identifyies repeats, 
# estimates mappability and finds all of the short scaffolds. 
# The output files include: 								
# ok.bed (regions to analyze in angsd etc)		
# map_repeat_summary.txt (summary of ref quality)							
#
####usage####
#User will need to input (paste) information for the following variables below:
#
#Genus-species: this is used in the directory naming as Erangi suggested, to make browsing
#a bit more doable for us humans
#accession: this is also used in the directory, to keep multiple reference assemblies separate as Black suggested
#n: the number of cpus you allocated in the SBATCH command above
#
#Example of defined variables below 
#
#genus_species=Balaenoptera-musculus
#accession=GCF_009873245.2
#n=32
#Once these have been defined, save and close slurrm job and submit
#sbatch /scratch/bell/$USER/theta/step1_download.sh
#
####end usage and notes####


###########################################
#ENTER INFORMATION FOR FOLLOWING VARIABLES#
###########################################

genus_species=
accession=
n=

########################
#DO NOT EDIT BELOW CODE#
########################

#move to species/accession directory

cd /scratch/bell/dewoody/theta/${genus_species}/

####assess mappability of reference####
genmap index -F ${accession}_ref/ref_100kb.fa -I index -S 50 # build an index 

# compute mappability, k = kmer of 100bp, E = # two mismatches
mkdir mappability
genmap map -K 100 -E 2 -T $n -I index -O mappability -t -w -bg                

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


# make list with the >100kb scaffolds -> duplicate with the step below making "ref_100kb.info"
#awk '{ print $1, $2, $2 }' ${accession}_ref/ref_100kb.fa.fai > ${accession}_ref/chrs.info

# replace column 2 with zeros -> duplicate with the step below making "ref_100kb.bed"
#awk '$2="0"' ${accession}_ref/chrs.info > ${accession}_ref/chrs.bed

# make tab delimited -> duplicate with the step below making "ref_100kb.bed"
#sed -i 's/ /\t/g' ${accession}_ref/chrs.bed

# make chrs.txt -> moved to below
#cut -f 1 ${accession}_ref/chrs.bed > ${accession}_ref/chrs.txt

# identify and remove sex-linked scaffolds with SATC -> remove SATC step from here
#make directory to house SATC files
#mkdir ${accession}_satc
#cd ${accession}_satc

# make bamfilelist
#ls ../sra/final_bams/*bam > bamlist
#sed -i 's/.bam//g' bamlist

# make a dir for the SATC output
#mkdir idxstats

# for each bam file calculate idxstats which reports alignment summary statistics
#ls -1 ../sra/final_bams/*bam | sed 's/\//\'$'\t/g' | cut -f 4| sed 's/_sorted.bam//g' > bamlist
#ls -1 ../sra/final_bams/*bam | sed 's/\//\'$'\t/g' | cut -f 4| sed 's/_sorted.bam//g' | while read -r LINE
#do
# Uncomment below if there are not any *bai files in directory
#samtools index ../sra/final_bams/${LINE}.bam
# idxstats
#samtools idxstats ../sra/final_bams/${LINE} > ${LINE}.idxstats
#cp ${LINE}.idxstats idxstats/${LINE}.idxstats
#done

# move bamlist to idxstats DIR and make list of idxstats files
#cp bamlist idxstats/bamlist
#cd idxstats
#ls *.idxstats > idxstats.txt

# run Sex Assignment Through Coverage (SATC)
#SATC="/scratch/bell/dewoody/SATC/satc.R"

# provide output prefix
#OUT="satc"

# list of idxstats files
#IN="idxstats.txt"

# setting minimum length to 100kb as we are not going to use scaffolds shorter than that

# choose between default
#Rscript --vanilla $SATC -i $IN -o $OUT --minLength 100000

# or 'useMedian'
#Rscript --vanilla $SATC -i $IN -o $OUT --useMedian TRUE --minLength 100000
# SATC will give a warning if 'useMedian' should be used

# SATC produce a PCA with individual sex assignment and boxplot with sex-linked scaffolds "satc_PCA_and_boxplot.png"
# combine all scaffolds both the ZX and the dodgy ones
#cat satc_sexlinked_scaff.list satc_XZ_scaff.list | sort -t _ -k 2 -g | uniq -u > sex.scafs 

#cd .. 

#sort ../${accession}_ref/chrs.txt > ../${accession}_ref/sorted_chrs.txt

# remove sex chromosome scaffolds from scaffold list
#comm -1 -3 ./idxstats/sex.scafs ../${accession}_ref/sorted_chrs.txt > ./autosomes.txt

#xargs samtools faidx ../${accession}_ref/ref_100kb.fa < ./autosomes.txt > ../${accession}_ref/autosomes_100kb.fa 

#samtools faidx ../${accession}_ref/autosomes_100kb.fa 

#move to species directory
#cd /scratch/bell/dewoody/theta/${genus_species}/ -> remove SATC step to here

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
cd /scratch/bell/${USER}/theta/source/
Rscript qc_reference_stats.R --args /scratch/bell/dewoody/theta/${genus_species}/ ${genus_species} ${accession} 
cd /scratch/bell/dewoody/theta/${genus_species}/

map=$(sed -n '1p' okmap.txt)
norepeats=$(sed -n '1p' norepeat.txt)
okbed=$(sed -n '1p' okbed.txt)

echo -e "${genus_species}\t $accession\t $map\t $norepeats\t $okbed" >> map_repeat_summary.txt

# make list of the bamfiles and index each file in theta directory for step5 in advance
# if you remove some individuals depending on their mapping rate, depth, or breadth, edit this file accordingly before step5
mkdir /scratch/bell/dewoody/theta/${genus_species}/theta/
cd /scratch/bell/dewoody/theta/${genus_species}/theta/
ls /scratch/bell/dewoody/theta/${genus_species}/sra/final_bams/*.bam > ./bam.filelist

# END
