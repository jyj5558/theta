#!/bin/bash
ml bioinfo
ml sra-toolkit
echo "٩(^‿^)۶"
echo Hello, what is the name of the accession ofr your species genome?
read species
echo Moving to directory that will house $species SRA files for FOUR channel chemistry
cd ./$species/sra/four_ch

echo " - - - -" ; echo " - - - -" ; echo " - - - -"         
echo Please specify a comma seperated list in a single row containing SRA for $species using FOUR channel chemistry
read SRA_4
echo $SRA_4 |  tr "," "\n" | while read g

do

echo downloading ${g}
echo ${g} | xargs prefetch --max-size 500GB ${g} -O ./fetch
echo ./fetch/${g}.sra | xargs fasterq-dump -e 6 -t ./ -O ./raw

done
echo "SRA four-channel chemistry files for $species have been downloaded and split. Please check to ensure all downloaded correctly"
