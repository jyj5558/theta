#!/bin/bash

#load modules
ml bioinfo BBMap/37.93
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo "٩(^‿^)۶"
echo Hello, what is the name of the accesion for the species genome you are about to download?
read species
echo Creating directory structure for $species
mkdir -p ./$species/${species}_ref
mkdir ./$species/${species}_rm
mkdir ./$species/${species}_gtf

mkdir ./$species/sra
mkdir ./$species/sra/two_ch
mkdir ./$species/sra/four_ch

mkdir ./$species/sra/two_ch/raw
mkdir ./$species/sra/two_ch/cleaned
mkdir ./$species/sra/two_ch/aligned

mkdir ./$species/sra/four_ch/raw
mkdir ./$species/sra/four_ch/cleaned
mkdir ./$species/sra/four_ch/aligned

echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo Please specify the path to the FTP refseq genome:
read ref

echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo Please specify the path to the FTP refseq rm out file:
read rm

echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo Please specify the path to the FTP refseq annotations in GTF format:
read gtf    

echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo Downloading reference genome, repeat masker out and annotation files for $species
wget $ref -O ./$species/${species}_ref/${species}.genomic.fna.gz -o ./$species/${species}_ref/log_ref 
wget $rm -O ./$species/${species}_rm/${species}.rm.out.gz -o ./$species/${species}_rm/log_rm
wget $gtf -O ./$species/${species}_gtf/${species}.gtf.gz -o ./$species/${species}_gtf/log_gtf

echo Downloads are now complete, decompressing files
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
gzip -d $species/*/*.gz
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2

echo Please confirm that all downloaded files were fully downloaded 100% below
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
cat ./$species/${species}_*/log* | grep "100%" 
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2

echo Please double check that all downloaded files are named correctly below and are not empty
ls -lh $species/GCF*
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo Now, the reference will be searched for a mitochondrion sequence
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
grep "mitochondrion" ./$species/${species}_ref/${species}.genomic.fna | cut -f 1 
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo If a sequence was printed out above, this will now be removed from the reference genome
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
filterbyname.sh include=f in=./$species/${species}_ref/${species}.genomic.fna out=./$species/${species}_ref/tmp.fasta names="mitochondrion" ow=t substring=t
rm ./$species/${species}_ref/${species}.genomic.fna
mv ./$species/${species}_ref/tmp.fasta ./$species/${species}_ref/${species}.genomic.fna
echo " - - - -" ;sleep 2; echo " - - - -" ;sleep 2 ; echo " - - - -" ;sleep 2
echo DONE-please check all files before moving on to the SRR download step
echo " "
echo "┏(-_-)┛┗(-_-<feff> )┓┗(-_-)┛┏(-_-)┓"
