#!/bin/bash
#SBATCH --job-name=TAR_genus-species
#SBATCH -A standby
#SBATCH -t 4:00:00 
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH -o %x_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=

#########################################################################
# this script archive, compress, and save Theta data files into Fortress				#
# the following file will be saved in Fortress:	  	#
# ${genus_species}.tar.gz 												#
#s cript was written to be run with SLURM job scheduler
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
cd /path/to/theta/${genus_species}/
rm ${genus_species}.tar.gz
echo "tarball started"
tar -cvf ${genus_species}.tar sra/final_bams/ ok.bed theta/ *ref/ref.fa *ref/ref.fa.fai *ref/chrs.txt *rm/repeats.bed
pigz ${genus_species}.tar -p 128
echo "tarball finished"

# Transfer to Fortress
echo "Transferring to Fortress started"
hsi -v put ${genus_species}.tar.gz /group/fnrdewoody/theta/${genus_species}.tar.gz
echo "Transferring to Fortress finished"

# Check if the file exists in the target Fortress directory "with the same memory" as in the local cluster directory 
echo "File size in Fortress is: "
echo `hsi -q ls -l /group/fnrdewoody/theta/${genus_species}.tar.gz`
echo "File size in Scratch is: "
echo `ls -l ./${genus_species}.tar.gz`

# And then remove original files by uncommenting below
#cd /path/to/theta/${genus_species}/
#rm ${genus_species}.tar.gz
#rm -r theta/
#rm -r sra/
#rm ok.bed
