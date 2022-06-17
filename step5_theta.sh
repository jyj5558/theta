#!/bin/bash
#SBATCH --job-name=S5_Genus-species
#SBATCH -A fnrquail
#SBATCH -t 10-00:00:00 
#SBATCH -N 1 
#SBATCH -n 64 
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="your email address (e.g.jeon96@purdue.edu) without quotation marks"


module --force purge
module load biocontainers
module load angsd
module load r
module load bcftools
module load htslib

#########################################################################
# this script identifies / filteres for usable sites and estimates theta#
# the following files to use in downstream analyses are created:	#
#########################################################################

####notes and usage####
#
##Notes##
# add target species "genus-species"
#Example:
#genus_species=Marmota-marmota-marmota
#usage:
#/scratch/bell/$USER/theta/step5_theta.sh
#
#
####end usage and notes####

##########################
#designate target species#
##########################

genus_species=

########Do not edit beneath this row########

#Path to parent directory of genus-species
PD=/scratch/bell/dewoody/theta/${genus_species}/

#First move to parent directory to be able to set variable MIND
cd $PD

##Designate min number of individuals, set to total numb of downloaded individuals divided by two
MIND=$((`wc -l < $PD/${genus_species}_SRA.txt` / 2))

#Just some absolute  paths to shorten commands
FINAL=/scratch/bell/dewoody/theta/${genus_species}/sra/final_bams
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta

#Move to directory which will house angsd/theta files

mkdir $THETA
cd $THETA
 
# make list of the bamfiles and index each file
ls ${FINAL}/*.bam > ./bam.filelist

# convert bed file to angsd format
awk '{print $1"\t"$2+1"\t"$3}' $PD/ok.bed > ./angsd.file

# index file
angsd sites index ./angsd.file

# estimate GL. Increase -P if you can/should employ more CPUs.
echo "Genotype likelihood estimation started"
angsd -P 64 -bam ./bam.filelist -sites ./angsd.file -anc $PD/*_ref/ref.fa \
-ref $PD/*_ref/ref.fa -dosaf 1 -gl 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 30 -minQ 30 -rf $PD/*_ref/chrs.txt -minInd $MIND -out out
echo "Genotype likelihood estimation done"

# obtain ML estimate of SFS using the folded realSFS. Increase -P if you can/should employ more CPUs.
echo "SFS estimation started"
realSFS -P 64 out.saf.idx  -fold 1 > out.sfs
echo "SFS estimation done"

# calculate theta for each site. Increase -P if you can/should employ more CPUs.
realSFS saf2theta  -P 64 out.saf.idx -sfs out.sfs -outname out

# estimate 
thetaStat print out.thetas.idx > out.thetas_persite.txt
echo "per-site Theta file created"

#Sliding window estimate
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
echo "window-based Theta file created"

# column 4 has Wattersons, column 9 has Tajima's D, and column 10 has Fu & Li's Fs
awk '{print $4}' theta.thetasWindow.gz.pestPG > Watterson
awk '{print $9}' theta.thetasWindow.gz.pestPG > TajimaD
awk '{print $10}' theta.thetasWindow.gz.pestPG > FuF

# get mean
meanW=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' Watterson)
meanD=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' TajimaD)
meanF=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' FuF)

# get SD
sdW=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanW2 += delta * ($1 - avg); } END { print sqrt(meanW2 / NR); }' Watterson)
sdD=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanD2 += delta * ($1 - avg); } END { print sqrt(meanD2 / NR); }' TajimaD)
sdF=$(awk '{delta = $1 - avg; avg += delta / NR; \
meanF2 += delta * ($1 - avg); } END { print sqrt(meanF2 / NR); }' FuF)

# print to file
echo -e "$PWD\t $meanW\t $sdW" \
>> Wattersons_theta_${genus_species}.txt
echo -e "$PWD\t $meanD\t $sdD" \
>> TajimaD_${genus_species}.txt
echo -e "$PWD\t $meanF\t $sdF" \
>> FuF_${genus_species}.txt
echo "Theta, D, F txt created"

####################################
# heterozygosity for each individual
####################################
#Increase "-P"s if you can/should employ more CPUs.

mkdir ./HET
OUTDIR='HET'

ls -1 ../sra/final_bams/*bam | sed 's/\//\'$'\t/g' | cut -f 4| sed 's/.bam//g' | while read -r LINE

do
echo "${LINE} heterozygosity estimation started"

angsd -i ../sra/final_bams/${LINE}.bam -ref $PD/*_ref/ref.fa -anc $PD/*_ref/ref.fa -dosaf 1 -rf $PD/*_ref/chrs.txt -sites ./angsd.file \
-minMapQ 30 -minQ 30 -P 64 -out ${OUTDIR}/${LINE} -only_proper_pairs 1 -baq 2 \
-GL 2 -doMajorMinor 1 -doCounts 1 -setMinDepthInd 5 -uniqueOnly 1 -remove_bads 1 

realSFS -P 64 -fold 1 ${OUTDIR}/${LINE}.saf.idx > ${OUTDIR}/est.ml

echo "${LINE} heterozygosity is" >> ./Het
cd $OUTDIR
Rscript -e 'a<-scan("est.ml"); a[2]/sum(a)' >>  ../Het
mv est.ml ${LINE}_est.ml
echo "${LINE} heterozygosity estimation done"

cd ../
done

####################################
# heterozygosity for population
####################################

# get mean
meanH=$(awk 'BEGIN{s=0;}{s=s+$2;}END{print s/NR;}' Het)

# get SD
sdH=$(awk '{delta = $2 - avg; avg += delta / NR; \
meanH2 += delta * ($2 - avg); } END { print sqrt(meanH2 / NR); }' Het)

# print to file
echo -e "$PWD\t $meanH\t $sdH" \
>> Het_${genus_species}.txt
echo "Population heterozygosity txt created"

#######
# ROHs
#######

angsd -b ./bam.filelist -dobcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -P 40

echo "ROH estimation started"
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' angsdput.bcf | bgzip -c > ${genus_species}.freqs.tab.gz
tabix -s1 -b2 -e2 ${genus_species}.freqs.tab.gz
bcftools roh --AF-file ${genus_species}.freqs.tab.gz --output ROH_${genus_species}.txt --threads 40 angsdput.bcf
echo "ROH txt created"

# END
