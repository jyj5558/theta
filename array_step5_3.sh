#!/bin/bash
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

###############################################################################
# this script identifies / filteres for usable sites and estimates theta, etc.#
# the following result files are created: 
# Het_${genus_species}.txt
###############################################################################

########################
#DO NOT EDIT BELOW CODE#
########################

# define variables
genus_species=$1

# move to directory which will house files
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta
cd $THETA

# calculate heterozygosity for population

# get mean
meanH=$(awk 'BEGIN{s=0;}{s=s+$2;}END{print s/NR;}' Het)

# get SD
sdH=$(awk '{delta = $2 - avg; avg += delta / NR; \
meanH2 += delta * ($2 - avg); } END { print sqrt(meanH2 / NR); }' Het)

# print to file
echo -e "$PWD\t $meanH\t $sdH" \
>> Het_${genus_species}.txt
echo "Population heterozygosity .txt created"

# END
