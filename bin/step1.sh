#!/bin/bash

module purge
module load biocontainers
module load samtools
module load bbmap

####notes and usage####
#
# This script downloads reference, and repeat 
# if a masked genome isn't available (i.e. rm.out), script will create one using the mammal database
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
    
#Move to the scratch space
cd $CLUSTER_SCRATCH/theta

####create directories and download reference genome, repeat masker, and annotation####
mkdir -p ./${genus_species}/${accession}_ref
mkdir -p ./${genus_species}/${accession}_rm
mkdir -p ./${genus_species}/${accession}_gtf
cd ${genus_species}

#Download reference genome
cd ${accession}_ref
wget ${pathway}${accession}_${assembly}_genomic.fna.gz 
gunzip ${accession}_${assembly}_genomic.fna.gz
cp ${accession}_${assembly}_genomic.fna original.fa # keep a copy of the original reference
cd ../

#Download repeatmasker file (if available)
cd ${accession}_rm
wget ${pathway}${accession}_${assembly}_rm.out.gz 
gunzip ${accession}_${assembly}_rm.out.gz
cp ${accession}_${assembly}_rm.out rm.out # keep a copy of the original repeatmasker
cd ../

#Print out file sizes for checking later
ls -lh ${accession}* > download_log

#Search for and remove mito seq in ref genome. Will only work if marked in assembly!
grep "mitochondrion" ${accession}_ref/original.fa | cut -f 1 > mito_header.txt #If no mitochondrial sequence present in reference, will be blank. Otherwise contain header
filterbyname.sh include=f in=${accession}_ref/original.fa out=${accession}_ref/original.tmp.fa names="mitochondrion" ow=t substring=t
rm ${accession}_ref/original.fa
mv ${accession}_ref/original.tmp.fa ${accession}_ref/original.fa      

###prep reference genome for mapping####
#Reduce fasta header length
reformat.sh in=${accession}_ref/original.fa out=${accession}_ref/new.fa trd=t -Xmx20g overwrite=T
#Sort by length
sortbyname.sh in=${accession}_ref/new.fa out=${accession}_ref/ref.fa -Xmx20g length descending overwrite=T
#Remove sequences smaller that 100kb prior to any repeatmasking
reformat.sh minlength=100000 in=${accession}_ref/ref.fa out=${accession}_ref/ref_100kb.fa overwrite=T
rm ${accession}_ref/new.fa

#Index ref.fa and ref_100kb.fa for step3, step4, and step5 
samtools faidx ${accession}_ref/ref_100kb.fa
samtools faidx ${accession}_ref/ref.fa

#Prep repeatmasked file for later processing, create a rm.out if one is not available. 
#Move into rm directory
cd ${accession}_rm/ 

#FILE1=$"rm.out"
#if [ -s $FILE1 ]; then
#Parse rm.out out to create bed coordinates
#  cat rm.out |tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' > repeats.bed # make bed file
#else	
#If no rm.out file is available, run RepeatMasker. Note, samtools conflict so had to purge first
#  module --force purge
# module load biocontainers
#  module load repeatmasker
#  RepeatMasker -pa 32 -a -qq -species mammals -dir . ../${accession}_ref/ref_100kb.fa 
#  cat ref_100kb.fa.out  | tail -n +4 | awk '{print $5,$6,$7,$11}' | sed 's/ /\t/g' > repeats.bed 
#fi

#Move back to species/accession directory
cd ../ 

#END