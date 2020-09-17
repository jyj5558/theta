# Theta project

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife.

August 2020

Author: Anna Brüniche-Olsen

## Steps for analyzing each species

Follow the pipeline below to qc reference assembly, qc read files, map reads, analyses etc.

QC of reference assembly
- download genome from NCBI (Black we need your scripts here)
- qc_reference.sh

QC and mapping of SRAs
- download SRAs from EBI-EMBL (Black we need your scripts here)
- trimgalore.sh
- mapping.sh
-realignment_single.sh (for individuals with a single SRA)
-realignment_multiple.sh (for individuals with multiple SRAs)

Analyses
- angsd.sh



