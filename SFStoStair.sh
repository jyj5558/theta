#!/bin/bash
#SBATCH --job-name=stair_genus-species
#SBATCH -A fnrpupfish
#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

module --force purge
module load biocontainers
module load angsd

genus_species=

PD=/scratch/bell/dewoody/theta/${genus_species}
THETA=/scratch/bell/dewoody/theta/${genus_species}/theta

cd $THETA

# designate min number of individuals, set to total numb of bam-filed individuals divided by two
MIND=$((`wc -l < ./bam.filelist` / 2))

# convert bed file to angsd format
awk '{print $1"\t"$2+1"\t"$3}' $PD/ok.bed > ./angsd.file

# index file
angsd sites index ./angsd.file

# estimate SFS; or just use out.sfs generated at step 5
echo "Site frequency spectrum estimation started"
angsd -bam ./bam.filelist -ref $PD/*_ref/ref.fa -anc $PD/*_ref/ref.fa -rf $PD/*_ref/chrs.txt -sites ./angsd.file \
-dosaf 1 -GL 2 \
-minInd $MIND -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 \
-out stair -P 64 

realSFS -P 64 stair.saf.idx  -fold 1 > stair.sfs
echo "Site frequency spectrum estimation finished"

cat stair.sfs | tr ' ' '\n' | awk -v OFMT=%.17g '{sum+=$1} END {print sum} # total number of observed nucleic sites, including polymorphic and monomorphic (needed for stairway plot)
# manually copy and paste numbers in ".sfa" file to ".blueprint" file after making a new ".blueprint" file, e.g., "Orcinus_orca_fold.blueprint"/ 
# For 'n' diploid samples, the site frequency spectrum (SFS) is the (2n+1) vector containing the proportion of site carrying 'k'-mutations. This means that the first element in the SFS is the proportion of sites where we don't observe any mutations, The second value is the proportion of sites where we observe 1 mutations. The last value is the proportion of sites we only observe mutations. It follows that the first and last column are the invariable categories. For folded sfs, the last half values will be "0".
# For "nrand", we suggest to use 4 numbers, roughly (nseq-2)/4, (nseq-2)/2, (nseq-2)*3/4, nseq-2. The range of the numbers is between 1 and nseq-2. 
# refer to Stairplot2 manual for detailed description of blueprint file.
cd /depot/fnrdewoody/apps/stairway_plot_v2.1.1
#cp two-epoch_fold.blueprint ${genus_species}_fold.blueprint
java -cp stairway_plot_es Stairbuilder ${genus_species}_fold.blueprint
bash ${genus_species}_fold.blueprint.sh
