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
- realignment for individuals with a single or with multiple SRAs
	-realignment_single.sh
	-realignment_multiple.sh

Analyses
- angsd.sh



