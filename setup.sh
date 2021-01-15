#!/bin/bash
echo Hello, what is the name of the species genome you are about to download?
read species
echo Creating directory structure for $species, please confirm that it was created correctly

mkdir -p ./$species/${species}_ref
mkdir ./$species/${species}_rm
mkdir ./$species/${species}_gtf
mkdir ./$species/sra
mkdir ./$species/sra/raw
mkdir ./$species/sra/cleaned
mkdir ./$species/sra/aligned
