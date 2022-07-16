#!/bin/bash
#SBATCH --job-name=TAR_genus-species
#SBATCH -A highmem
#SBATCH -t 1-00:00:00 
#SBATCH -N 1 
#SBATCH -n 128
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=

#########################################################################
# this script archive, compress, and save Theta data files into Fortress				#
# the following file will be saved in Fortress:	  	#
# ${genus_species}.tar.gz 												#
#########################################################################
####notes and usage####
#
##Notes##
# add target species "genus-species"
#Example:
#genus_species=Marmota-marmota-marmota
#usage:
#/scratch/bell/$USER/theta/Theta_tarball.sh
#
#
####end usage and notes####

##########################
#designate target species#
##########################

genus_species=

########################
#Do not edit code below#
########################

# Tarball
cd /scratch/bell/dewoody/theta/${genus_species}/
echo "tarball started"
tar -cvf ${genus_species}.tar sra/final_bams/ ok.bed theta/ 
pigz ${genus_species}.tar -p 128
echo "tarball finished"

# Transfer to Fortress
echo "Transferring to Fortress started"
hsi -v put ${genus_species}.tar.gz /group/fnrdewoody/theta/${genus_species}.tar.gz
echo "Transferring to Fortress finished.

# Manually check if the file exists in the target Fortress directory "with the same (or larger) memory" as in the local cluster directory by uncommenting below on your terminal screen (not by job submission):
#hsi -q ls -l "/group/fnrdewoody/theta/${genus_species}.tar.gz"

# And then remove original files by uncommenting below
#cd /scratch/bell/dewoody/theta/${genus_species}/
#rm ${genus_species}.tar.gz
#rm -r theta/
#rm -r sra/
#rm ok.bed
