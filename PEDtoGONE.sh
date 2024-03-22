#!/bin/bash
#SBATCH --job-name=GONE_genus-species
#SBATCH -A 
#SBATCH -t 14-00:00:00 
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

#script was written to be run with SLURM job scheduler

module --force purge
module load biocontainers
module load bioinfo
module load angsd
module load plink/1.90b6.4

genus_species=
n=

########Do not edit beneath this row########

# absolute paths of parent and theta directory of genus-species shorten commands
PD=/path/to/theta/${genus_species}
THETA=/path/to/theta/${genus_species}/theta

cd $THETA

# designate min number of individuals, set to total numb of bam-filed individuals divided by two
MIND=$((`wc -l < ./bam.filelist` / 2))
#IND=$((`wc -l < ./bam.filelist`))
#LN=$((`wc -l < $PD/*_ref/chrs.txt`/10))
#awk -v ln=$LN 'NR <= ln' $PD/*_ref/chrs.txt > chrs-plink.txt # subset 10% longest contigs for GONE
#awk -v ln=$LN 'NR <= 1' $PD/*_ref/chrs.txt > chrs-plink2.txt # subset the longest contigs for GONE
awk 'NR <= 200 {print $1}' $PD/*_ref/ref_100kb.fa.fai > chrs-plink.txt # ref_100kb.fa.fai can be used instead of chrs.txt

# plink file for GONE
angsd -bam bam.filelist -ref $PD/*_ref/ref.fa -rf chrs-plink.txt -sites ./angsd.file -out plink_GONE \
-doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 2 -doCounts 1 -doMaf 1 -postCutoff 0.99 -SNP_pval 1e-6 -geno_minDepth 5 \
-minMapQ 30 -minQ 30 -minInd $MIND -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -baq 2 -setMinDepthInd 5 -P $n

# convert to plink .ped and .map format
plink --tfile plink_GONE --recode --allow-extra-chr --out ${genus_species}_GONE

# create a chromosome code file
awk 'BEGIN {OFS=" "}; {print $0, NR}' chrs-plink.txt > chr-codes.txt

# save a copy of .map file
cp ${genus_species}_GONE.map ${genus_species}_GONEcp.map

# change chromosome code to numerical values
awk -vOFS=" " 'FNR==NR{arr[$1]=$2;next}{$1=arr[$1]; print}' chr-codes.txt ${genus_species}_GONE.map > ${genus_species}_GONE2.map
mv ${genus_species}_GONE2.map ${genus_species}_GONE.map
#join -o '2.2 1.2 1.3 1.4' $genus_species.map chr-codes.txt

# copy to the directory having GONE script
cp ${genus_species}_GONE.ped /path/to/GONE-master/GONE-master/Linux/${genus_species}_GONE.ped 
cp ${genus_species}_GONE.map /path/to/GONE-master/GONE-master/Linux/${genus_species}_GONE.map 

# run GONE
cd /dpath/to/GONE-master/GONE-master/Linux/
bash script_GONE.sh ${genus_species}_GONE

# move back result files to species directory
mv *_${genus_species}* $THETA

cd $THETA

# designate min number of individuals, set to total numb of bam-filed individuals divided by two
MIND=$((`wc -l < ./bam.filelist` / 2))
IND=$((`wc -l < ./bam.filelist`))
LN=$((`wc -l < $PD/*_ref/chrs.txt`/10))
#awk -v ln=$LN 'NR <= ln' $PD/*_ref/chrs.txt > chrs-plink.txt # subset 10% longest contigs for GONE
awk -v ln=$LN 'NR <= 1' $PD/*_ref/chrs.txt > chrs-plink2.txt # subset the longest contigs for GONE

# plink file for GONE
angsd -bam bam.filelist -ref $PD/*_ref/ref.fa -rf chrs-plink2.txt -sites ./angsd.file -out plink_GONE2 \
-doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 2 -doCounts 1 -doMaf 1 -postCutoff 0.99 -SNP_pval 1e-6 -geno_minDepth 5 \
-minMapQ 30 -minQ 30 -minInd $IND -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -baq 2 -P $n

# convert to plink .ped and .map format
plink --tfile plink_GONE2 --recode --allow-extra-chr --out ${genus_species}_GONE2

# create a chromosome code file
awk 'BEGIN {OFS=" "}; {print $0, NR}' chrs-plink2.txt > chr-codes2.txt

# save a copy of .map file
cp ${genus_species}_GONE2.map ${genus_species}_GONE2cp.map

# change chromosome code to numerical values
awk -vOFS=" " 'FNR==NR{arr[$1]=$2;next}{$1=arr[$1]; print}' chr-codes2.txt ${genus_species}_GONE2.map > ${genus_species}_GONE3.map
mv ${genus_species}_GONE3.map ${genus_species}_GONE2.map
#join -o '2.2 1.2 1.3 1.4' $genus_species.map chr-codes.txt

# copy to the directory having GONE script
cp $genus_species.ped /path/to/GONE-master/GONE-master/Linux/${genus_species}_GONE2.ped 
cp $genus_species.map /path/to/GONE-master/GONE-master/Linux/${genus_species}_GONE2.map 

# run GONE
cd /path/to/GONE-master/GONE-master/Linux/
bash script_GONE.sh ${genus_species}_GONE2

# move back result files to species directory
mv *_${genus_species}* $THETA

cd $THETA

# END
