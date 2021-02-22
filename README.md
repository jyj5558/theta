# Theta project - Feb 2021

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife populations. NOTE, these scripts assume that a SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).

## Project organization links
[Google doc with data and links](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469); 
[Weekly agenda and notes](https://docs.google.com/document/d/1vyvKtTTdbAaev23nXTlfw-awJjivq9ENdmS1YWzZW4I/) 

## Script and directory structure overview
Script overview and order
Step 1: Reference download and QC (step1_ref_download_QC.sh) 
Step 2: Population data identification, download, and cleaning (step2_SRA_download_clean.sh)
Step 3: ...
- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)

QC of bam dataset
- qc_bams.sh & quantile_thresholds.R

Analyses
- angsd.sh

Directory overview
To make things easier to automate later, we will use a standardized directory structure based on the target species scientific name (Genus-species), genome assembly accession number (accession), and other labels. The structure looks like this:

Genus-species

-accession_ref: reference genome assembly

-accession_rm: repeat-masked reference genome assembly

-accession_gtf: reference genome assembly annotation

-mappability: reference genome assembly quality assessment

-sra: holds several levels of population data

--raw: raw SRAs, kept for now to help with debugging

--cleaned: cleaned SRAs that are ready for mapping

--mapped: SRAs mapped to reference genome

## Reference download and QC - step1_ref_download_QC.sh
This script downloads reference, repeat, and annotation data and then identifyies repeats, estimates mappability and finds all of the short scaffolds. The output files include: 	
1.ref.fa (reference file with scaffolds>100kb)							
2.ok.bed (regions to analyze in angsd etc)		
3.map_repeat_summary.txt (summary of ref quality)							
4.if a masked genome isn't available (i.e. rm.out), script will create one using the mammal repeat library --- we should update this if we move on from mammals!

Usage: step1_download_QC.sh Genus-species accession PUid pathway
Genus-species: this is used in the directory naming as Erangi suggested, to make browsing  a bit more doable for us humans
accession: this refers to the NCBI ref seq accession number and is also used in the directory, to keep multiple reference assemblies separate as Black suggested
PUid: this tells this script where to find the R scripts (is there a better way to do this?)
pathway: include NCBI path up to but not including file extension, e.g. /genomes/all/GCF/001/890/085/GCF_001890085.1_ASM189008v1/GCF_001890085.1_ASM189008v1

eventually, we will write a script to automate running these but for now we will test one species at a time. To run this step, use something like this:
sbatch /scratch/bell/jwillou/theta/step1_ref_download_QC.sh Manis-javanica GCF_001685135.1 /genomes/all/GCF/001/890/085/GCF_001890085.1_ASM189008v1/GCF_001890085.1_ASM189008v1 jwillou

Note that this requires genmap. I've installed this in DeWoody scratch and then added a PATH update for this location, but we should watch for errors on this step.


## Population data identification, download, and cleaning
###Identifying population data
We need to find sequence data for populations corresponding to the reference genome assemblies that are available. To do this, we will first look for SRA data by species and then check to see that it meets the criteria outlined below.

1.Start looking for population data at NCBI (https://www.ncbi.nlm.nih.gov/genbank/), EMBL (https://www.ebi.ac.uk/), publications (Google scholar etc.), or in Andrew Black’s downloaded SRA directory (available in this repo under SRA_metadata) for the target species. The target species is listed in [our Google doc of data](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469). Be sure to check the species carefully, sometimes the target species is the host, meaning that they will show up in the search but did not actually get sequenced.
Andrew Black's [SRA_metadata](./SRA_metadata/) contains separate files for many of the species. You can find BioProject accession number, BioSample number ..etc from the respective file to use and then double check those numbers in NCBI.

2.Once you find a sequencing record for the target species, check that it is whole-genome resequencing data from Illumina technologies Hi-Seq or Nova Seq but not MiSeq. Record the platform used in the google doc so we can look for systematic biases later. This is usually found in the NCBI record but if not you will need to look in the reference noted in the NCBI record.

3.Be sure that the read length should is 150 base pairs or longer and that the reads should are paired-end. (Paired end reads will have R1 and R2 files.)

4.Check that the sequenced individuals were captured from wild populations. (We will exclude laboratory-bred or other captive-bred populations.) 

5.Most of the papers have a section named ‘data availability’ and this section has the BioProject accession number(PRJNAxx…) which you can follow to find SRAs. (The SRAs are the unique identifiers for the NCBI-submitted data we want.) Follow the Supplementary information section of the paper for more information if it is available.

6.Per species, SRAs from a maximum of 8 individuals should be collected. The 8 individuals should belong to one population (as defined by the original authors) If there are multiple populations use sequences from the population with 8 individuals. If all the populations are less than 8 individuals, use data from the population with the highest number of individuals sequenced.
For an instance, if there are 2 populations of black bears, one with 6 and the other with 4 individuals, record SRAs for the population with 6 individuals. We need at least 2 individuals per species (minimum number) to include in the project.

6.Make sure each SRA belongs to one individual and 8 SRAs belong to 8 different individuals. There can be multiple SRAs for one individual and these multiple SRAs can be included in the project. The pipeline is developed with a space to add multiple SRAs. 

### Download and clean population data - step2_SRA_download_clean.sh
This script downloads SRAs for our target species and then checks the initial quality of the reads, cleans the reads, and checks the read quality again after cleaning. The output files include:
1.SRAs downloaded and unzipped, in sras/raw
2.SRAs cleaned (trimmed), in sras/cleaned

Usage: step2_SRA_download_clean.sh Genus-species 
Genus-species: This is used in the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.
SRAs.txt: A comma separated list of SRAs to download (should be on a single row). Note that SRAS.txt must be in /scratch/bell/dewoody/theta/${genus_species}/.



