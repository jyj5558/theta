# Theta project - since Feb 2021 (for reference to future script users), published at _Evolutionary Applications_ Aug 2024

These scripts are developed for investigating genomic diversity and its association with threat categories in wildlife populations. 
NOTE, these scripts assume that Purdue's SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).
When running jobs, if the jobs stay still for a longer time than you expected, try again with a larger memory or larger number of CPUs (e.g., "#SBATCH --mem=" or #SBATCH -n" arguments).

There is a Nextflow pipeline script which you can submit as a slurm job in a bash script (or just "nextflow run main_theta.nf") with an associated input csv file example.
Note that all script files should be housed under "/bin/" directory within the directory that houses the nextflow script. 
I recommend housing it in "<your scratch>/theta" directory since the nextflow will generate many and large intermediate files.

Below is explaining what each step does (step1 to step5) when you run individually.

## Directory structure overview

To make things easier to automate, we will use a standardized directory structure based on the target species scientific name (Genus-species), genome assembly accession number (accession), and other labels. The folder structure looks like this:

      -Genus-species/

             -accession_ref: reference genome assembly

             -accession_rm: repeat-masked reference genome assembly

             -accession_gtf: reference genome assembly annotation

             -mappability: reference genome assembly quality assessment
             
             -sra: raw, cleaned, aligned sra files, and ready-to-analyze bam files
             
             -theta: ANGSD output GD metrics

All of the Genus-species directories will be created at <your scratch>/theta

## Cloning Git repo to your cluster directory
Clone the theta git repos to user scratch path (e.g., /scratch/bell/blackan/) by running the following code
```
cd <your scratch>
git clone https://github.com/jyj5558/theta
```


## Step 1. Reference downloading
**Identifying reference assembly data**

We need to find reference assembly accession information for as many species as possible. Whenever possible, we will target Chromosome level RefSeq assemblies on NCBI. Once the target assembly is located for a particular species, we will use the accession number, assembly name and url path (found under the FTP links on the right side of NCBI) and assembly name to download resources effeciently.

**Download reference assembly - step1.sh**

This script downloads reference, repeat, and annotation data.  	

If a masked genome isn't available (i.e. rm.out), script will create one using the mammal repeat library --- you should update this if you use non-mammals!

User will need to input (paste) information for the following variables within the step1 file (the one in the archive folder) if using a stand-alone script (i.e., not using nextflow script):

`Genus_species`: this is used in the directory naming to make browsing a bit more doable for us humans

`accession`: this is also used in the directory, to keep multiple reference assemblies separate

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
#SBATCH --job-name=S1_Balaenoptera-musculus

Once these have been defined, save and close slurm job script

Before running, make sure that the program genmap is installed at $USER $HOME directory and in $PATH, such as:
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

Submit slurm job (when running a stand-alone script)
```
sbatch <your scratch>/theta/step1_ref_download.sh
```

>**Quality check #1:**
>Inside the target species directory, confirm that repeat masking (or extraction of rm.out file) worked before proceeding to step2!
```
head *_rm/repeats.bed
```

## Step 2. Population data identification, download, and cleaning
**Setting up for SRA download**

Before the user can download SRA files for the first time, they need to make sure that "remote access is enabled". This should only need to be done once:

Load modules:
```
ml biocontainers
ml sra-tools
```
And now enable access
```
vdb-config -i graphical
```
A panel will pop up and enter "E" to enable, "O" to click okay, and "S" to save and close out of screen. 
The two tools we will be using (prefetch and fasterq-dump) will download to the working directory and remove temporary files afterwards. If a bottleneck is encountered (space wise) we may need to add a new path to a directory which will house large temporary files.

**Identifying population data**

We (the authors of paper) found sequence data for a single population corresponding to the reference genome assemblies that are available. For each species with an existing reference assembly, we looked for short read sequencing data that has been generated for a single population. In particular, the data sets need to have these characteristics:

1. whole-genome resequencing data (not RNAseq, RADseq, etc., also check if it is exom sequencing data.)
2. generated using the Illlumina Genome Analyzer, Hi-Seq, NovaSeq, and NextSeq platforms (not MiSeq) (Library should not be constructed for Reduced Representation or anything that substantially lowers "Bases")
3. be comprised of paired-end reads with lengths of ideally 150 base pairs or longer (paired end reads will  have R1 and R2 files), but we tried to download regardless of base pair lengths and record the read length info in SRA metadata file.
4. sequenced individuals were captured from wild populations (we excluded laboratory-bred or other captive-bred populations). We did not include non-wild populations. 
5. sequenced individuals should all belong to the same population; this can be tough to identify for some species, so we went with the definition used by the authors that published the sequencing data (assuming there is an associated publication)
6. a minimum of 2 individuals but ideally up to 25 individuals were sequenced independently (more than one individual cannot be contained in a single sequencing project); if there are  multiple files (SRAs) for one individual, chose the best one based on "Bases".


**Download and clean population data - step2.sh**

This script downloads SRAs for our target species and then checks the initial quality of the reads, cleans the reads, and checks the read quality again after cleaning. The output files include:
1. SRAs downloaded and unzipped, in sra/raw
2. SRAs cleaned (trimmed), in sra/cleaned


Usage:

User will need to enter the following information within the step2 file (the one in the archive folder) if using a stand-alone script (i.e., not using nextflow script):
`Genus-species`: This will need to be enterred in the the above script and follow the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.

`NOTE` make sure that the Genus-species-SRAs.txt is updated first: Make sure you have the species SRA metadata in the correct format and in the theta directory BEFORE running. It should be in theta/SRA_metadata/genus-species.txt

`NOTE`, add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
```
#SBATCH --job-name=S2_Balaenoptera-musculus
```

Now, when ready to run, simply submit the following command:
```
sbatch <your scratch>/theta/step2_SRA_download_clean.sh
```
>**Quality check #2:**
>Inside the target species directory, confirm that you have the right number of downloaded raw / cleaned samples. There should be the same number of >samples as the number of rows (-1) in the Species_sra file. 
Inside the target species directory, run the following code:
```
cat *SRA.txt | tail -n +2 | wc -l
```
The numerical output should match that obtained from the following:
```
ls -1 sra/cleaned/*1.fq | wc -l
```
If the numbers do not match, figure out which sample(s) did not download, why and download individually if necessary. If the number match, proceed to next quality check

>**Quality check #3:**
>Inside the target species directory, print out the statistics from trimgalore. If any samples have less than 80% of reads being retained, note this sample (to remove from downstream analysis)

```
grep -r "Reads written (passing filters):" sra/cleaned/*report.txt
```
After noted, proceed to step 3.

## Steps 3: Mapping cleaned SRA fastq files to the reference assembly
**- step3.sh** 

This script will align each paired fastq file to the reference, validate the output sam file, sort based upon read coordinate, mark pcr duplicates and perform indel realignment. Mapping and depth statistics will be output. 
-Binary alignment mapping (.bam) files: Curated files will be found within ./final_bams
-Summary statistic files will be found within: ./alignment

Usage:

User will need to enter the following information within the step3 file (the one in the archive folder) if using a stand-alone script (i.e., not using nextflow script):
`Genus-species`: This will need to be enterred in the the above script and follow the directory naming (see step 1), so we need it to make sure the SRAs are saved in the correct location.

`NOTE` make sure that the Genus-species-SRAs.txt is updated first: Make sure you have the species SRA metadata in the correct format and in the theta directory BEFORE running. It should be in theta/SRA_metadata/genus-species.txt

`NOTE`, add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
```
#SBATCH --job-name=S3_Balaenoptera-musculus
```

When ready to run, simply submit the following command:
```
sbatch <your scratch>/theta/step3_mapping.sh
```

>**Quality check #4:**
>Inside the target species directory, evaluate the mapping rates. If any outliers are identified (e.g., low mapping rate (i.e. less than 80%), low depth (i.e. less than 1x), low breadth (i.e., less than 80%)) note these for potential removal during the last step. Also, confirm that there are the correct number of rows (below) as there are samples. A discrepency could identify a problem sample. 
 
 Alignment rates:
 ```
 grep -r "+ 0 mapped (" sra/final_bams/*mapping.txt
 ```
Depth:
```
cat sra/final_bams/*depth*
```
Breadth:
```
cat sra/final_bams/*breadth*
```
Record individuals' depth and breadth in ascending order of SRR numbers.
Again, note any problem samples. Re-align if necessary before proceeding to step4:
To realign, the easiest way is to modify the script and change out the 'for loop' for the sample ID (e.g. CHANGE "for i in `ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq`" TO "for i in SRR10596315 SRR10596317" if you want to realign just SRR10596315 and SRR10596317), then resave script as a new name and run.
If the number of samples excluding the problematic samples is larger than 2, then just move on to the next step with notes on the samples.

## Step 4: QC the reference. Create ok.bed file, containing mappable sites, non-repeat sites, and autosomal sites
**- step4_qc_noSATC.sh**

This script identifyies repeats, estimates mappability and finds all of the short scaffolds. 

The output files include: 			

1. ok.bed (regions to analyze in angsd etc)		
2. map_repeat_summary.txt (summary of ref quality)
3. chr.txt (chromosome names to limit regions and facilitate the analysis in angsd etc)		


User will need to enter the following information within the step4 file (the one in the archive folder) if using a stand-alone script (i.e., not using nextflow script):
`Genus-species`: 
`accession`: 

`NOTE`, add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
```
#SBATCH --job-name=S4_Balaenoptera-musculus
```
Once entered, simply submit the following command:
```
sbatch <your scratch>/theta/step4_qc_noSATC.sh
```
>**Quality check #5:**
>Inside the target species directory, Confirm the final bed file is produced (and not empty) and that the sex scaffolds script worked correctly.

Bed file should have three columns (scaffold, start and stop)
```
head ok.bed
```

## Step 5: Estimate thetas, ROH, and individual heterozygosity from the alignment files
**- step5_theta_v2.sh**

This script will estimate the Site Frequency Spectrum and calculate Watterson's Theta (in sliding windows), Nucleotide diveristy, Tajima's D, Fu & Li's Fs, individual heterozygosities, population heterozygosity, and runs of homozygosity. As ANGSD will take a lot of memory, you might need to reserve the whole half node (64 CPUs; or even the whole node - 128 cores). The output file will contain mean with SD (except ROH).
The output file:
1. Thetas_${Genus-species}.txt
2. Ind_Het_${Genus-species}.txt (containing individual heterozygosity values; 'comma'-separated)
3. Pop_Het_${Genus-species}.txt (containing the population-level heterozygosity value)
4. ROH_${Genus-species}_PL.txt (ROH values based on phred-score scaled genotype likelihood)

If some samples were dropped due to quality issues from upstream analyses, $MIND in step5 script should be changed accordingly (by editing the bam.filelist after making it manually or renaming the problematic sample.bam in advance to something not to be included in bam.filelist (e.g., sample.nobam) are the easiest).

Usage:

User will need to enter the following information within the step5 file (the one in the archive folder) if using a stand-alone script (i.e., not using nextflow script):
`Genus-species`:
`accession`:

`NOTE`, add the name of the genus-species to the job name as well (for debugging/tracking purposes) as in:
```
#SBATCH --job-name=S5_Balaenoptera-musculus
```

To run, simply submit the following command:
```
sbatch <your scratch>/theta/step5_theta_v2.sh
```

Record population-level Nucleotide diversity, Watterson's theta, Tajima's D, Fu and Li's F, Heterozygosity, F(ROH) >= 100kb, and F(ROH) >= 1mb.
Also record individual-level Heterozygosity in ascending order of SRR numbers (as the same order of individual depth and breadth for statistical correlation test; the samples should have been automatiaclly sorted and processed in the ascending order in step3 and so in downstream steps, so nothing needs to be done to sort them, but just double-check.).

#################################
#################################


