# Theta project

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife.

August 2020

## Steps for analyzing each species

Follow the pipeline below to qc reference assembly, qc read files, map reads, analyse data etc.

QC of reference assembly
- qc_reference.sh & repeatmasker_names.R
- qc_reference_stats.sh & qc_reference_stats.R

QC and mapping of SRAs
- trimgalore.sh
- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)

#### pipeline is tested until this step ####
QC of bam dataset
- qc_bams.sh & quantile_thresholds.R



Analyses
- angsd.sh

test edit


