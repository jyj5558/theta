#!/bin/bash
#SBATCH --job-name=mapSites
#SBATCH -A standby
#SBATCH -t 4:00:00 
#SBATCH -N 1 
#SBATCH -n 1 

module load bioinfo
module load r

# cd $SLURM_SUBMIT_DIR

cat species.txt | while read -r LINE

do

cd $LINE

Rscript /scratch/snyder/a/abruenic/scripts/qc_reference_stats.R

map=$(sed -n '1p' okmap.txt)
repeat=$(sed -n '1p' norepeat.txt)
okbed=$(sed -n '1p' okbed.txt)

echo -e "$LINE\t $map\t $repeat\t $okbed" \
>> /scratch/snyder/a/abruenic/pupfish/map_repeat_summary.txt

cd ..

done

#END
