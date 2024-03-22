#!/bin/bash
#SBATCH --job-name=S1_Genus-species
#SBATCH -A 
#SBATCH -t 12-00:00:00 
#SBATCH -N 1 
#SBATCH -n 32
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="your email address (e.g.jeon96@purdue.edu) without quotation marks"

module load bioinfo
module load bioawk
module load samtools
module load BBMap

####notes####
#
# This script downloads reference, repeat, and annotation data
#if a masked genome isn't available (i.e. rm.out), script will create one using the mammal 
#scripts was written to be run with SLURM job scheduler
#
####usage####
#
#User will need to input (paste) information for the following variables below:
#
#Genus-species: this is used in the directory naming as Erangi suggested, to make browsing
#a bit more doable for us humans
#accession: this is also used in the directory, to keep multiple reference assemblies
#separate as Black suggested
#pathway: include full NCBI url to FTP site (containing directory)                                          
#assembly: name of assembly
#n: the number of cpus you allocated in the SBATCH command above
#
#Example of defined variables below 
#genus_species=Balaenoptera-musculus
#accession=GCF_009873245.2
#pathway=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/
#assembly=mBalMus1.pri.v3
#n=32
#Once these have been defined, save and close slurrm job and submit
#sbatch /scratch/bell/$USER/theta/step1_download.sh
####end usage and notes####

###########################################
#ENTER INFORMATION FOR FOLLOWING VARIABLES#
###########################################

genus_species=
accession=
pathway=
assembly=
n=

########################
#DO NOT EDIT BELOW CODE#
########################

#Move to the scratch space
cd $CLUSTER_SCRATCH

####create directories and download reference genome, repeat masker, and annotation####
mkdir -p ./$genus_species/${accession}_ref
mkdir ./$genus_species/${accession}_rm
mkdir ./$genus_species/${accession}_gtf
cd $genus_species

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

#Download annotation file (if available)
cd ${accession}_gtf
wget ${pathway}${accession}_${assembly}_genomic.gtf.gz 
gunzip ${accession}_${assembly}_genomic.gtf.gz
cp ${accession}_${assembly}_genomic.gtf gtf.gtf # keep a copy of the original annotation
cd ../

#print out file sizes for checking later
ls -lh ${accession}* > download_log

#search for and remove mito seq in ref genome. Will only work if marked in assembly!
grep "mitochondrion" ${accession}_ref/original.fa | cut -f 1 > mito_header.txt #If no mitochondrial sequence present in reference, will be blank. Otherwise contain header
filterbyname.sh include=f in=${accession}_ref/original.fa out=${accession}_ref/original.tmp.fa names="mitochondrion" ow=t substring=t
rm ${accession}_ref/original.fa
mv ${accession}_ref/original.tmp.fa ${accession}_ref/original.fa      

###prep reference genome for mapping####
#Reduce fasta header length
reformat.sh in=${accession}_ref/original.fa out=${accession}_ref/new.fa trd=t -Xmx20g overwrite=T
#sort by length
sortbyname.sh in=${accession}_ref/new.fa out=${accession}_ref/ref.fa -Xmx20g length descending overwrite=T
#remove sequences smaller that 100kb prior to any repeatmasking
bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' ${accession}_ref/ref.fa > ${accession}_ref/ref_100kb.fa
rm ${accession}_ref/new.fa

#index ref.fa and ref_100kb.fa for step3, step4, and step5 
samtools faidx ${accession}_ref/ref_100kb.fa
samtools faidx ${accession}_ref/ref.fa

#prep repeatmasked file for later processing, create a rm.out if one is not available. 
#move into rm directory
cd ${accession}_rm/ 

FILE1=$"rm.out"
if [ -s $FILE1 ]
then
	# parse rm.out out to create bed coordinates
	cat rm.out |tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' > repeats.bed # make bed file
else	
	# if no rm.out file is available run RepeatMasker. Note, samtools conflict so had to purge first
	module --force purge
	module load biocontainers/default
	module load repeatmasker
	RepeatMasker -pa $n -a -qq -species mammals -dir . ../${accession}_ref/ref_100kb.fa 
	cat ref_100kb.fa.out  | tail -n +4 | awk '{print $5,$6,$7,$11}' | sed 's/ /\t/g' > repeats.bed 
fi

#move back to species/accession directory
cd ../ 

####prep annotation file for later processing####
cd ${accession}_gtf/ #move into gtf directory

FILE1=$"gtf.gtf"
if [ -s $FILE1 ]
then
	# Print message
	echo "gtf file properlly formatted" > log_gtf
else	
	# Print message
	printf "no annotation available" > log_gtf
fi
#DONE

