# Theta project - Feb 2021

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife populations. NOTE, these scripts assume that a SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).

## Project organization links
- [Google doc with data and links](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469)
- [Weekly agenda and notes](https://docs.google.com/document/d/1vyvKtTTdbAaev23nXTlfw-awJjivq9ENdmS1YWzZW4I/) 

## Directory structure overview

To make things easier to automate later, we will use a standardized directory structure based on the target species scientific name (Genus-species), genome assembly accession number (accession), and other labels. The folder structure looks like this:

Genus-species/
     - accession_ref: reference genome assembly
     - accession_rm: repeat-masked reference genome assembly
     - accession_gtf: reference genome assembly annotation
     - mappability: reference genome assembly quality assessment

All of the Genus-species directories will be created at /scratch/bell/dewoody/theta

## Cloning Git repo to your cluster directory
Locate to your $CLUSTER_SCRATCH (e.g., /scratch/bell/blackan/) and run the code below to clone Theta git repo 
```
git clone https://github.com/AnnaBrunicheOlsen/theta
```

## Step 1. Reference downloading and QC
**Identifying reference assembly data**

We need to find reference assembly accession information for as many species as possible. Whenever possible, we will target Chromosome level RefSeq assemblies on NCBI. Once the target assembly is located for a particular species, we will use the accession number and file path (found under the FTP links on the right side of NCBI) and assembly name to download resources effeciently.

**Download and perform QC on reference assembly - step1_ref_download_QC.sh**
This script downloads reference, repeat, and annotation data and then identifyies repeats, estimates mappability and finds all of the short scaffolds. The output files include: 	
1. ref.fa (reference file with scaffolds>100kb)							
2. ok.bed (regions to analyze in angsd etc)		
3. map_repeat_summary.txt (summary of ref quality)							
4. if a masked genome isn't available (i.e. rm.out), script will create one using the mammal repeat library --- we should update this if we move on from mammals!

User will need to input (paste) information for the following variables within the _step1_ref_download_QC.sh file:_

`Genus-species`: this is used in the directory naming as Erangi suggested, to make browsing
a bit more doable for us humans

`accession`: this is also used in the directory, to keep multiple reference assemblies
separate as Black suggested

`pathway`: include full NCBI url to FTP site (containing directory)

`assembly`: name of genome assembly

Example of defined variables below
```
genus_species=Balaenoptera-musculus
accession=GCF_009873245.2
pathway=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/
assembly=mBalMus1.pri.v3
```
`NOTE`, it is strongly recommended to add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
#SBATCH --job-name=qcREF-Balaenoptera-musculus

Once these have been defined, save and close slurm job script

Before running, make sure that the program genmap is installed at $USER home directory and in $PATH, such as:
```
export PATH=$PATH:~/genmap-build/bin
```
```
genmap --help
```
If program usage prints out, you are ready to submit the SLURMM job (below). If the command is not found, set $PATH accordingly.
To install genmap in your home directory, (refer to https://github.com/cpockrandt/genmap for details)
```
$ cd ~
$ git clone --recursive https://github.com/cpockrandt/genmap.git
$ mkdir genmap-build && cd genmap-build
$ module load cmake
$ cmake ../genmap -DCMAKE_BUILD_TYPE=Release
$ make genmap
$ ./bin/genmap
```

Submit slurm job
```
sbatch /scratch/bell/$USER/theta/step1_ref_download_QC.sh
```

Eventually, we will write a script to automate running these in batch but for now we will download one species at a time.


## Step 2. Population data identification, download, and cleaning
**Identifying population data**

We need to find sequence data for populations corresponding to the reference genome assemblies that are available. For each species with an existing reference assembly (see [our Google doc of data](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469) for more information), we are looking for short read sequencing data that has been generated for a single population. In particular, the data sets need to have these characteristics:

1. whole-genome resequencing data (not RNAseq, RADseq, etc., also check if it is exom sequencing data.)
2. generated using the Illlumina Hi-Seq or Nova Seq platforms (not MiSeq) (Library should not be constructed for Reduced Representation or anything that substantially lowers "Bases")
3. be comprised of paired=end reads with lengths of 150 base pairs or longer (paired end reads will  have R1 and R2 files)
4. sequenced individuals were captured from wild populations (we will exclude laboratory-bred or other captive-bred populations; but please record this info, too, in "Notes")
5. sequenced individuals should all belong to the same population; this can be tough to identify for some species, so we will go with the definition used by the authors that published the sequencing data
6. a minimum of 8 individuals but ideally up to 25 individuals were sequenced independently (more than one individual cannot be contained in a single sequencing project); it is fine if there are  multiple files (SRAs) for one individual

There are two pathways to find these data, it is important to use both methods for all target species.

1) 
a: Open up the [Google doc sheet](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469), locate and copy the name of your species of interest 

b: Paste this into the [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi), remove the "-" and search. Copy the search result as a hyperlink into the [Google doc sheet](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469) under genus-species.

c: Click on the "SRA Experiments" number, which will list the number of Entrez records associated with this species. Click on all SRA experiments of relevance (i.e., Illumina WGS paired end data) and click the "Send to" button, choosing clipboard. Repeat for each relevant BioProject until clipboard is populated with all available usable data. (When there are too many SRAs, apply "advanced search filter" right below the search window to narrow down data. Add relevant filters, such as "Add 'paired'" & "Add 'WGS'". Click "Search" and at the search result page, click "Send results to Run selector" to send all SRAs at once.)

d: Navigate to clipboard and 'Send to' to the 'Run selector'. Click on the 'metadata' icon to download file (you can filter available data once more in this page). Alternatively, go up to 'Send to' at the top of the page, choose 'File' and 'Accessions List', then click 'Create File'

e: Open this file with Excel, convert text to columns (comma seperated delimiter) and remove/edit the file until it has the same order and headers as below:
(ReadLength should be included, and it can be found in SRA search page). It should be in tab seperated format

```
Run	BioProject	BioSample	ScientificName	pop	Instrument	LibraryLayout	sex AvgSpotlength Bases
ERR3299081	PRJEB8272	SAMEA5577702	Marmota marmota marmota	Gsies	NextSeq 500	PAIRED	female 300     5316224638
ERR3299082	PRJEB8272	SAMEA5577703	Marmota marmota marmota	Gsies	NextSeq 500	PAIRED	male 300  3316224671
ERR3299083	PRJEB8272	SAMEA5577704	Marmota marmota marmota	Gsies	NextSeq 500	PAIRED	male 300  3312124192
ERR3294906	PRJEB8272	SAMEA5577694	Marmota marmota marmota	Mauls	Illumina HiSeq 2500	PAIRED	female 300     4315224637
```

f: Copy and paste this into the species file in [SRA_metadata](./SRA_metadata/).

g: Download / pull git repository to $CLUSTER_SCRATCH or paste new metadata for your target species directly into existing Genus-species file in $USER theta folder

 
2) To verify that all available data was found for a given species, use google scholar to look for literature directly. For each genus-species, use search terms  "Whole Genome" OR "resequencing" OR "genomic" and >2010. Manually curate search results. Most of the papers have a section named ‘data availability’ and this section has the BioProject accession number(PRJNAxx…) which you can follow to find the data we need. Otherwise, on NCBI, search the project number written in SRA metadata file (e.g. PRJNA54005) targetting "All Databases", and you can find "Literature" box showing the number of associated publications in PubMed, PubMed Central, etc. Here you can find the project description, too, which normally has a statement about sample sources.

Fill out all missing fields for target species in the [Google doc](https://docs.google.com/spreadsheets/d/1u9Zxzcms1DdeV0k8qyJpFboO81r1Uvl8udIt8PRjUSk/edit#gid=235995469): 

`NOTE`: It is important to estimate coverage for resequencing data. We really want to target low coverage WGS. Do a back of the envelope calculation for each SRA file to ensure that it is less than 20x. Black found it easiest to do in the excel file when creating the metadata by) using the number of bases in the SRA/2 / genome size and averaging it among all SRA samples (just look for outliers above 20x). This information (i.e., average) can be put in the 'coverage' column in the shared google doc.

**Download and clean population data - step2_SRA_download_clean.sh**

This script downloads SRAs for our target species and then checks the initial quality of the reads, cleans the reads, and checks the read quality again after cleaning. The output files include:
1. SRAs downloaded and unzipped, in sra/raw
2. SRAs cleaned (trimmed), in sra/cleaned


Usage:

User will need to enter the following information:
`Genus-species`: This will need to be enterred in the the above script and follow the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.

`NOTE` make sure that the Genus-species-SRAs.txt is updated first: Make sure you have the species SRA metadata in the correct format and in the theta directory BEFORE running. It should be in theta/SRA_metadata/genus-species.txt

`NOTE`, it is strongly recommended to add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
#SBATCH --job-name=SRA-Balaenoptera-musculus

To run, simply submit the following command:
```
sbatch /scratch/bell/$USER/theta/step1_download_QC.sh
```

**Steps 3: Mapping cleaned SRA fastq files to the reference assembly - step3_mapping.sh** 

This script will align each paired fastq file to the reference, validate the output sam file, sort based upon read coordinate, mark pcr duplicates and perform indel realignment. Mapping and depth statistics will be output. 
-Binary alignment mapping files: Curated files will be found within ./final_bams
-Summary statistic files will be found within: ./alignment

Usage:

User will need to enter the following information:
`Genus-species`: This will need to be enterred in the the above script and follow the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.

`NOTE` make sure that the Genus-species-SRAs.txt is updated first: Make sure you have the species SRA metadata in the correct format and in the theta directory BEFORE running. It should be in theta/SRA_metadata/genus-species.txt

`NOTE`, it is strongly recommended to add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
#SBATCH --job-name=SRA-Balaenoptera-musculus

To run, simply submit the following command:
```
sbatch /scratch/bell/$USER/theta/step1_download_QC.sh
```

`NOTE`, please check for any bugs with the step3 script, it is still new!


**Steps 4: QC bams? - NEEDS WORK** 

QC of bam dataset
- qc_bams.sh & quantile_thresholds.R

Analyses
- angsd.sh

