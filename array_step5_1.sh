#!/bin/bash
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

###############################################################################
# this script identifies / filteres for usable sites and estimates theta, etc.#
# the following result files are created: 
# NucleotideDiversity_${genus_spceis}.txt
# WattersonsTheta_${genus_spceis}.txt
# TajimaD_${genus_species}.txt
# FuF_${genus_species}.txt
# ROH_${genus_species}.txt
# and many other interim files (e.g., out.thetas_persite.txt, angsdput.bcf)
###############################################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# load modules
module --force purge
module load biocontainers
module load angsd
module load bcftools
module load htslib

# define variables
genus_species=$1
n=$2
accession=$3

# path to parent directory of genus-species
PD=/scratch/bell/dewoody/theta/${genus_species}/

# first move to parent directory to be able to set variable MIND
cd $PD

# just some absolute  paths to shorten commands
FINAL=/scratch/bell/dewoody/theta/${genus_species}/sra/final_bams
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta

# move to directory which will house angsd/theta files
mkdir -p $THETA/HET 
cd $THETA
 
# make list of the bamfiles and index each file
# if you remove some individuals depending on their mapping rate, depth, or breadth, edit this file accordingly
ls $FINAL/*.bam > ./bam.filelist

# designate min number of individuals, set to total numb of bam-filed individuals divided by two
MIND=$((`wc -l < ./bam.filelist` / 2))

# convert bed file to angsd format
awk '{print $1"\t"$2+1"\t"$3}' $PD/ok.bed > ./angsd.file

# index file
angsd sites index ./angsd.file

# estimate GL
echo "Genotype likelihood estimation started"
angsd -P ${n} -bam ./bam.filelist -sites ./angsd.file -anc $PD/*_ref/ref.fa \
-ref $PD/*_ref/ref.fa -dosaf 1 -gl 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 30 -minQ 30 -rf $PD/*_ref/chrs.txt -minInd $MIND -out out
echo "Genotype likelihood estimation done"

# obtain ML estimate of SFS using the folded realSFS
echo "SFS estimation started"
realSFS -P ${n} out.saf.idx  -fold 1 > out.sfs
echo "SFS estimation done"

# calculate theta for each site
echo "Theta estimation started"
realSFS saf2theta -P ${n} out.saf.idx -sfs out.sfs -outname out

# estimate 
thetaStat print out.thetas.idx > out.thetas_persite.txt
echo "per-site Theta file created"

# sliding window estimate
thetaStat do_stat out.thetas.idx -win 50000 -step 10000 -outnames theta.thetasWindow.gz
echo "window-based Theta file created"

# column 4 has Wattersons, column 5 has Nucleotide diversity, column 9 has Tajima's D, and column 10 has Fu & Li's Fs
awk '{print $4}' theta.thetasWindow.gz.pestPG > Watterson
awk '{print $5}' theta.thetasWindow.gz.pestPG > NucDiv
awk '{print $9}' theta.thetasWindow.gz.pestPG > TajimaD
awk '{print $10}' theta.thetasWindow.gz.pestPG > FuF

# get mean
meanW=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' Watterson)
meanN=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' NucDiv)
meanD=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' TajimaD)
meanF=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' FuF)

# get SD
sdW=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanW2 += delta * ($1 - avg); } END { print sqrt(meanW2 / NR); }' Watterson)
sdN=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanN2 += delta * ($1 - avg); } END { print sqrt(meanN2 / NR); }' NucDiv)
sdD=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanD2 += delta * ($1 - avg); } END { print sqrt(meanD2 / NR); }' TajimaD)
sdF=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanF2 += delta * ($1 - avg); } END { print sqrt(meanF2 / NR); }' FuF)

# print to file
echo -e "$PWD\t $meanW\t $sdW" \
>> WattersonsTheta_${genus_species}.txt
echo -e "$PWD\t $meanN\t $sdN" \
>> NucleotideDiversity_${genus_species}.txt
echo -e "$PWD\t $meanD\t $sdD" \
>> TajimaD_${genus_species}.txt
echo -e "$PWD\t $meanF\t $sdF" \
>> FuF_${genus_species}.txt
echo "Pi, Theta, D, F .txt created"

#######
# ROHs
#######

echo "File conversion to bcf started"
angsd -b ./bam.filelist -dobcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -P ${n}
echo "bcf file created"

echo "ROH estimation started"
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' angsdput.bcf | bgzip -c > ${genus_species}.freqs.tab.gz
tabix -s1 -b2 -e2 ${genus_species}.freqs.tab.gz
bcftools roh --AF-file ${genus_species}.freqs.tab.gz --output ROH_${genus_species}_raw.txt --threads ${n} angsdput.bcf
echo "ROH raw file created"

echo "ROH raw file parsing started"
python $CLUSTER_SCRATCH/theta/ROHparser.py ${genus_species} ${accession}

# END
