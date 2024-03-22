#!/bin/bash
#SBATCH --job-name=S5_Genus-species
#SBATCH -A 
#SBATCH -t 14-00:00:00 
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module --force purge
module load biocontainers
module load angsd
module load r
module load bcftools
module load htslib

#########################################################################
# this script identifies / filteres for usable sites and estimates theta#
# the following result files are created:	
# Thetas_${genus_species}.txt, IndHet_${genus_species}.txt, PopHet_${genus_species}.txt,
# ROH_${genus_species}.txt
# script was written to be run with SLURM job scheduler
#########################################################################

####notes and usage####
#
##Notes##
#add target species "genus-species" and the number of cpus "n" you allocated in SBATCH command above
#Example:
#genus_species=Marmota-marmota-marmota
#accession=GCF_001458135.1
#n=64
#
#usage:
#/scratch/bell/$USER/theta/step5_theta_v2.sh
#
#
####end usage and notes####

##########################
#designate target species#
##########################

genus_species=
accession=
n=

########Do not edit beneath this row########

# absolute paths of parent and theta directory of genus-species shorten commands
PD=/path/to/theta/${genus_species}
THETA=/path/to/theta/${genus_species}/theta

# move to directory which will house angsd/theta files
mkdir -p $THETA
cd $THETA

# designate min number of individuals, set to total numb of bam-filed individuals divided by two
MIND=$((`wc -l < ./bam.filelist` / 2))

# convert bed file to angsd format
awk '{print $1"\t"$2+1"\t"$3}' $PD/ok.bed > ./angsd.file

# index file
angsd sites index ./angsd.file

# estimate GL
echo "Genotype likelihood estimation started"
angsd -bam ./bam.filelist -ref $PD/*_ref/ref.fa -anc $PD/*_ref/ref.fa -rf $PD/*_ref/chrs.txt -sites ./angsd.file \
-dosaf 1 -GL 2 -doMajorMinor 1 \
-minInd $MIND -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 \
-out out -P $n 
echo "Genotype likelihood estimation done"


# obtain ML estimate of SFS using the folded realSFS
echo "SFS estimation started"
realSFS -P $n out.saf.idx  -fold 1 > out.sfs
echo "SFS estimation done"

# calculate theta for each site
echo "Theta estimation started"
realSFS saf2theta out.saf.idx -sfs out.sfs -outname out -P $n 

# estimate 
thetaStat print out.thetas.idx > out.thetas_persite.txt
echo "per-site Theta file created"

#Sliding window estimate
thetaStat do_stat out.thetas.idx -win 50000 -step 50000 -outnames theta.thetasWin.gz
echo "window-based Theta file created"

# column 4 has Wattersons, column 5 has Nucleotide diversity, column 9 has Tajima's D, and column 10 has Fu & Li's Fs
awk '{print $4,$5,$9,$10,$14}' theta.thetasWin.gz.pestPG > Thetas

# get mean
sumSites=$(awk 'BEGIN{s=0;}{s=s+$5;}END{print s;}' Thetas)
sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas)
sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas)
meanW=$(awk "BEGIN {print $sumW/$sumSites}")
meanN=$(awk "BEGIN {print $sumN/$sumSites}")
meanD=$(awk 'BEGIN{d=0;}{d=d+$3;}END{print d/NR;}' Thetas)
meanF=$(awk 'BEGIN{f=0;}{f=f+$4;}END{print f/NR;}' Thetas)

echo -e "Genus_Species\tNucleotide_Diversity\tWatterson_Theta\tTajima_D\tFu&Li_F\n${genus_species}\t$meanN\t$meanW\t$meanD\t$meanF\t" >> Thetas_${genus_species}.txt
echo "Thetas (Pi, Theta, D, F) .txt created"

####################################
# heterozygosity for each individual #rerunning Loxodonta, Syncera, did not apply doCounts and setMinDepthInd filters
####################################

mkdir ./HET
OUTDIR='HET'

cat ./bam.filelist | sed 's/\//\'$'\t/g' | cut -f 9 | sed 's/.bam//g' | while read -r LINE

do
echo "${LINE} heterozygosity estimation started"

angsd -i ../sra/final_bams/${LINE}.bam -ref $PD/*_ref/ref.fa -anc $PD/*_ref/ref.fa -rf $PD/*_ref/chrs.txt -sites ./angsd.file \
-dosaf 1 -GL 2 -doMajorMinor 1 \
-minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -baq 2 -doCounts 1 -setMinDepthInd 5 \
-out ${OUTDIR}/${LINE} -P $n
 
realSFS -P $n -fold 1 ${OUTDIR}/${LINE}.saf.idx > ${OUTDIR}/${LINE}_est.ml

cd ${OUTDIR}
Rscript -e 'args<-commandArgs(TRUE); LINE<-args[1]; a<-scan(paste(LINE,"est.ml", sep="_")); a[2]/sum(a)' ${LINE} >>  ../Het
echo "${LINE} heterozygosity estimation done"
cd ../
done

cat ./Het | cut -d " " -f 2 | tr '\n' ',' > ./IndHet_${genus_species}.txt

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
>> PopHet_${genus_species}.txt
echo "Population heterozygosity .txt created"

#######
# ROHs
#######

echo "File conversion to bcf started"
angsd -b ./bam.filelist -ref $PD/*_ref/ref.fa -rf $PD/*_ref/chrs.txt -sites ./angsd.file \
-dobcf 1 -gl 2 -dopost 1 -domajorminor 1 -domaf 1 \
-minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -snp_pval 1e-6 -P $n
echo "bcf file created"

echo "ROH estimation started"
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' angsdput.bcf | bgzip -c > ${genus_species}.freqs.tab.gz
tabix -s1 -b2 -e2 ${genus_species}.freqs.tab.gz
bcftools roh --AF-file ${genus_species}.freqs.tab.gz --output ROH_${genus_species}_PLraw.txt --threads $n angsdput.bcf
echo "ROH raw files created"

echo "ROH raw file parsing started"
python3 $CLUSTER_SCRATCH/theta/ROHparser.py ${genus_species} ${accession}
# END
