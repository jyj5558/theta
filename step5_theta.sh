#!/bin/bash
#SBATCH --job-name=S4_Genus-species
#SBATCH -A fnrdewoody
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 40 
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module --force purge
module load biocontainers
module load angsd

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
#/scratch/bell/$USER/theta/step4_theta.sh
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

##Designate min number of individuals, set to total numb of downloaded individuals
MIND=`cat $PD/${genus_species}_SRA.txt | wc -l`

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

# estimate GL
angsd -P 40 -bam ./bam.filelist -sites ./angsd.file -anc $PD/*_ref/ref.fa -ref $PD/*_ref/ref.fa -dosaf 1 -gl 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -rf $PD/*_ref/chrs.txt -minInd 2 -out ./out

# obtain ML estimate of SFS using the folded realSFS
realSFS -P 40 out.saf.idx  -fold 1 > out.sfs

# calculate theta for each site
realSFS -P 40 saf2theta out.saf.idx -sfs out.sfs -outname out

# estimate 
#This is temporary until added to biocontainer path, but should work for now:
thetaStat do_stat out.thetas.idx

#This is temporary until added to biocontainer path, but should work for now:
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz

# column 4 has Wattersons
awk '{print $4}' theta.thetasWindow.gz > Watterson

# get mean
meanW=$(awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}' Watterson)

# get SD
sdW=$(awk '{delta = $1 - avg; avg += delta / NR; \
mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' Watterson)

# print to file
echo -e "$PWD\t $meanW\t $sdW" \
>> ${THETA}/Wattersons_theta_${Genus-species}.txt

# END
