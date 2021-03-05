# Theta project - Feb 2021

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife populations. NOTE, these scripts assume that a SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).

## Project organization links
- [Google doc with data and links](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469)
- [Weekly agenda and notes](https://docs.google.com/document/d/1vyvKtTTdbAaev23nXTlfw-awJjivq9ENdmS1YWzZW4I/) 

## Directory structure overview

To make things easier to automate later, we will use a standardized directory structure based on the target species scientific name (Genus-species), genome assembly accession number (accession), and other labels. The structure looks like this:

Genus-species
- accession_ref: reference genome assembly
- accession_rm: repeat-masked reference genome assembly
- accession_gtf: reference genome assembly annotation
- mappability: reference genome assembly quality assessment
- sra: holds several levels of population data
- sra/raw: raw SRAs, kept for now to help with debugging
- sra/cleaned: cleaned SRAs that are ready for mapping
- sra/mapped: SRAs mapped to reference genome

All of the Genus-species directories will be at /scratch/bell/dewoody/theta

## Step 1. Reference downloading and QC
**Identifying reference assembly data**

We need to find reference assembly accession information for as many species as possible. Whenever possible, we will target RefSeq assemblies on NCBI Genbank. Once the target assembly is located for a particular species, we will use the accession number and file path (found under the FTP links on the right side) to dowload assemblies effeciently.

**Download and perform QC on reference assembly - step1_ref_download_QC.sh**
This script downloads reference, repeat, and annotation data and then identifyies repeats, estimates mappability and finds all of the short scaffolds. The output files include: 	
1. ref.fa (reference file with scaffolds>100kb)							
2. ok.bed (regions to analyze in angsd etc)		
3. map_repeat_summary.txt (summary of ref quality)							
4. if a masked genome isn't available (i.e. rm.out), script will create one using the mammal repeat library --- we should update this if we move on from mammals!

Usage: step1_download_QC.sh Genus-species accession PUid pathway
- Genus-species: this is used in the directory naming as Erangi suggested, to make browsing  a bit more doable for us humans
- accession: this refers to the NCBI ref seq accession number and is also used in the directory, to keep multiple reference assemblies separate as Black suggested
- PUid: this tells this script where to find the R scripts (is there a better way to do this?)
- pathway: include NCBI path up to but not including file extension, e.g. /genomes/all/GCF/001/890/085/GCF_001890085.1_ASM189008v1/GCF_001890085.1_ASM189008v1

To run this step, use something like this:

sbatch /scratch/bell/jwillou/theta/step1_ref_download_QC.sh Manis-javanica GCF_001685135.1 /genomes/all/GCF/001/890/085/GCF_001890085.1_ASM189008v1/GCF_001890085.1_ASM189008v1 jwillou

Eventually, we will write a script to automate running these but for now we will test one species at a time. Note that this requires genmap; I've installed this in DeWoody scratch and then added a PATH update for this location, but we should watch for errors on this step.


## Step 2. Population data identification, download, and cleaning
**Identifying population data**

We need to find sequence data for populations corresponding to the reference genome assemblies that are available. For each species with an existing reference assembly (see [our Google doc of data](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469) for more information), we are looking for short read sequencing data that has been generated for a population. In particular, the data sets need to have these characteristics:

1. whole-genome resequencing data (not RNAseq, RADseq, etc.)
2. generated using the Illlumina Hi-Seq or Nova Seq platforms (not MiSeq) 
3. be comprised of paired=end reads with lengths of 150 base pairs or longer (paired end reads will  have R1 and R2 files)
4. sequenced individuals were captured from wild populations (we will exclude laboratory-bred or other captive-bred populations)
5. sequenced individuals should all belong to the same population; this can be tough to identify for some species, so we will go with the definition used by the authors that published the sequencing data
6. a minimum of 8 individuals but ideally up to 25 individuals were sequenced independently (more than one individual cannot be contained in a single sequencing data file/SRA); it is fine if there are  multiple files (SRAs) for one individual

There are several ways to find these data, but the workflow that seems to work the best often begins with NCBI and then transitions to primarily literature. (Other usable starting places include  EMBL (https://www.ebi.ac.uk/) and Andrew Black's [SRA_metadata](./SRA_metadata/), which contains information about NCBI files for many of the species of interest.)

1. Start by searching for your target species on NCBI (https://www.ncbi.nlm.nih.gov/genbank/), targeting BioProjects. B
2. Use some of the BioProject filters (e.g. genome sequencing, mammals, etc.) to limit the number of records for the particular target species. Then, read through the descriptions of the BioProjects to find one that seems to meet the criteria above. At this stage, it is important to notice what was sequenced (we want whole genomes and not RNA) and how many BioSamples are included in the project (we need at least 2 but 8 individuals is the target). Be sure to check the species of each BioProject carefully because sometimes our target species is the host of the sequenced species, meaning that they will show up in the search but will not have been sequenced and so will not be what we need. 
3. At this point it is often helpful to pull up the reference paper noted on the BioProjects page to see if the samples collected for this project will meet our criteria. If the samples will not work for us, go back to the list of BioProjects and try another one. 
4. Once you find a BioProject that looks like it meets the criteria we need, click through to the list of BioSamples (from the table near the bottom of the page) and look at the BioSamples that relate to that BioProject. Scan through the BioSamples until you find 8 (maximum) individuals that were sampled in the same area. Sometimes, sampling location is not noted on the BioSample page. Unless all of the samples for the project were collected from the same location (this information would be in the paper), missing sampling location information will mean that we cannot use these data and you'll need to start back at the BioProjects list again.
5. Once you have found 8 BioSamples that will work (or, if you found at least 2 and exhausted all other BioProject options), you can add the data we need to [our Google doc of data](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469). Note that if you cannot find what you need looking through NCBI, you can use google scholar to look for literature directly. Most of the papers have a section named ‘data availability’ and this section has the BioProject accession number(PRJNAxx…) which you can follow to find the data we need. 

Data needed in the [Google doc](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469): 

1. BioProject number (if relevant)
2. link to the reference paper noted in NCBI (or alternatively the one you used to find the data)
3. sequencing platform used; this is usually found in the NCBI record but if not you will need to look in the reference noted in the NCBI record
5. name, sampling location, or other identifying information about the population you chose so that we can look back later for more information as needed; it is simplest to use either the name used by the authors in the paper or the sampling location noted in the BioSample information
6. total number of populations available, if there is more than one set in the paper/BioProject
7. total number of individuals with data available in the population you chose
8. number of individuals with exported SRA or sample accession numbers added to the google sheet
9. comma separated list of SRAs or sample accession numbers; there should be no spaces between items in the list

Hint for extracting SRAs/BioSample accession numbers: In some cases, it will be much easier to extract sample accession numbers instead of all of the SRAs, and this is totally fine and will work in our pipeline. To do this effeciently, you can tick the boxes in the list of BioSamples for the sample numbers you want to export, go up to 'Send to' at the top of the page, choose 'File' and 'Accessions List', then click 'Create File'. This will generate a list of accession numbers in a text file, but they are newline separated instead of comma separated. Use find and replace to change this, and then paste the whole line into the google doc. See [this video](https://github.com/AnnaBrunicheOlsen/theta/blob/master/help/export%20sample%20numbers.mov) for a walkthrough.

**Download and clean population data - step2_SRA_download_clean.sh**

This script downloads SRAs for our target species and then checks the initial quality of the reads, cleans the reads, and checks the read quality again after cleaning. The output files include:
1. SRAs downloaded and unzipped, in sras/raw
2. SRAs cleaned (trimmed), in sras/cleaned

Usage: step2_SRA_download_clean.sh Genus-species 
- Genus-species: This is used in the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.
- SRAs.txt: A comma separated list of SRAs to download (should be on a single row). Note that SRAS.txt must be in /scratch/bell/dewoody/theta/${genus_species}/.


**Steps 3+:** To be completed.
- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)

QC of bam dataset
- qc_bams.sh & quantile_thresholds.R

Analyses
- angsd.sh

