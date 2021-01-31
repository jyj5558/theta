#!/bin/bash
ml sra-toolkit
echo "٩(^‿^)۶"
echo Hello, what is the name of the species genome?
read species
echo Moving to directory that will house $species SRA files
cd ./$species/sra/

echo " - - - -" ; echo " - - - -" ; echo " - - - -"         
echo Please specify a comma seperated list in a single row containing SRA for $species
read SRA
echo $SRA |  tr "," "\n" | while read i
do
echo downloading $i
echo ${i} | xargs prefetch --max-size 200GB $i -O ./fetch
echo ./fetch/${i}.sra | xargs fasterq-dump -e 6 -t ./ -O ./raw -b 10G -c 50G

done
echo "SRA files for $species have been downloaded and split. Please check to ensure all downloaded correctly"
