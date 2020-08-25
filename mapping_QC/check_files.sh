cat birds.txt  | while read -r LINE

do

cd $LINE

echo $LINE

file=$"diploid_raw.pdf"
if [ -s "$file" ]
then 
   echo " diploid_raw.pdf OK "
else
   echo " diploid_raw.pdf missing "
fi

file=$"diploid_repeatmasker.pdf"
if [ -s "$file" ]
then 
   echo " diploid_repeatmasker.pdf OK "
else
   echo " diploid_repeatmasker.pdf missing"
fi

cd ..

done


cat birds.txt  | while read -r LINE

do

cd $LINE

echo $LINE

file=$"SNPs.vcf"
if [ -s "$file" ]
then 
   echo " SNPs.vcf OK "
else
   echo " SNPs.vcf missing "
fi

cd ..

done


cat birds.txt  | while read -r LINE
do
cd $LINE
echo $LINE
size=$(grep "average length" slurm-52*.out)
echo $size
size=$(grep "insert size average" slurm-52*.out)
echo $size
cd ..
done