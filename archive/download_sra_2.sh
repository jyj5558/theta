#!/bin/bash
ml bioinfo
ml sra-toolkit
echo "٩(^‿^)۶"
echo Hello, what is the name of the accession ofr your species genome?
read species
echo Moving to directory that will house $species SRA files for TWO channel chemistry
cd ./$species/sra/two_ch

echo " - - - -" ; echo " - - - -" ; echo " - - - -"         
echo Please specify a comma seperated list in a single row containing SRA for $species using TWO channel chemistry
read SRA_2
echo $SRA_2 |  tr "," "\n" | while read i

do
echo downloading $i
echo ${i} | xargs prefetch --max-size 500GB $i -O ./fetch
echo ./fetch/${i}.sra | xargs fasterq-dump -e 6 -t ./ -O ./raw

done
echo "SRA two-channel chemistry files for $species have been downloaded and split. Please check to ensure all downloaded correctly"
