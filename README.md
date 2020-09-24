# Theta project

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife.

August 2020

Author: Anna Brüniche-Olsen

## Steps for analyzing each species

Follow the pipeline below to qc reference assembly, qc read files, map reads, analyses etc.

QC of reference assembly
- qc_reference.sh
- qc_reference_stats.sh & qc_reference_stats.R

QC and mapping of SRAs
- download_SRAs.sh
- trimgalore.sh
- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)

Identify sexlinked scaffolds, excessive heterozygosity and error rates estimation


Analyses
- angsd.sh



