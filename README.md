# Theta project
Jan 2021

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife populations. NOTE, these scripts assume that a SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).

## Steps for analyzing each species (one at a time)

STEP 1: Create folder structure and download reference genome, repeater masker out file and annotations (if available). NOTE, below the Cheetah genome is used as an example only. User's will need to designate the genome assembly for the species they are targetting. 

a) For a current list of all mammalian refseq assemblies, go to: https://www.ncbi.nlm.nih.gov/datasets/genomes/?txid=40674&source_db=RefSeq
and locate your species of interest. 

b) Identify the best / most current assembly (based upon N50, assembly levlel, etc) and click on the refseq hyperlink beginning with "GCF_", which will link to the NCBI genome summary web-page. 

Once here, locate and navigate to the "FTP directory for RefSeq assembly" in the right hand column. This will pull up the File Transfer Protocal web page where the reference genome (e.g. GCF_002007445.1_ASM200744v2_genomic.fna.gz), repeat masker out (e.g. GCF_002007445.1_ASM200744v2_rm.out.gz) and annotations (e.g. GCF_002007445.1_ASM200744v2_genomic.gtf.gz) are available for download. 

Copy the name of the reference genome build, and the link addresses to all three files (if available)

Now that you have this information recorded, run the [setup.sh](./setup.sh) script in the ./theta folder. The user will be prompted to enter the following information, which will create the necessary folder structure, download all files and decompress them. If the user species does not have one of these listed files, just enter NA when prompted).

Now, run the [setup.sh](./setup.sh) script in the ./theta folder to 1) auto-create the proper directory structure and b) download and decompress all files:

For example:

```
bash ./setup.py

Hello, what is the name of the species genome you are about to download?
GCF_002007445.1_ASM200744v2

Creating directory structure for GCF_002007445.1_ASM200744v2, please confirm that it was created correctly

Please specify the path to the FTP refseq genome:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_genomic.fna.gz

Please specify the path to the FTP refseq rm out file:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_rm.out.gz

Please specify the path to the FTP refseq annotations:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_genomic.gtf.gz

```

Once completed, check the log files to ensure that all fully downloaded (100%), decompressed, named correctly and are not empty

```
cat log* | grep "100%" 
#There should be one line printed out for each file downloaded

ls */
#All files should be decompressed
```

#STOP HERE, UNTESTED BELOW!!!!

STEP 2:

Now that the files have been download and decompressed, the reference needs to be QC'd using the following scripts:

  a) [qc_reference.sh](./qc_reference.sh) & [repeatmasker_names.R](./repeatmasker_names.R)
  
  b) [qc_reference_stats.sh](./qc_reference_stats.sh) & [qc_reference_stats.R](./qc_reference_stats.R)


STEP 3:

Locate and download SRA files for population level reads:

This document outlines the steps that can be taken to download Sequence Read Archive (SRA) files. The code listed below are demonstrated using my (Andrew Black’s) scratch space and paths. You will need to modify these accordingly. 
While you can download SRA files by submitting a job, which makes more sense if downloading a lot of files, we will be doing this logged into an interactive session:
```
sinteractive -A standby -t 4:00:00 -n 6
```

Now, we need to load two modules:
```
ml bioinfo sra-toolkit/2.10.0
```
The next step is to set the path to where you want to store the temporary (large) sequence read archive files. By default, these go into your home directory, and they can quickly use up all of your home directory space. Therefore, we want to specify a location to store them that is on your scratch drive (as you have terabytes of free space here!)
The location that you will want to specify (to store these SRA files) is dependant on the user path, but for example:
```
/scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/
```
Now that we have determined the path to where you will want to save these files, copy the path above. Once copied, run the following command: 
```
vdb-config -i --interactive-mode textual
```
This interactive session will ask you to specify some options. 
Type 4, and then paste the path that you previously copied above and hit <ENTER>.
After you input this information Type Y and <ENTER>.

Next, we’ll need to download the SRA accession(s) that you are interested in.

Copy the list of SRR accessions you are interested into a file called species_x_srr.txt 

Let’s put this file at the following path:
```
/scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/SRR.txt
```
This file should look something like this:
```
cat ./SRR.txt

SRR6656231
SRR6656187
SRR6656230
SRR6656186
SRR6656233
SRR6656185
. . .
```
Now that we have the text file created, we need to use the “prefetch command” which will download the SRR accessions to the computing infrastructure. 
Now, you should navigate to the folder that contains the species_x_srr.txt  file, and type this command:
```
cat ./SRR.txt | xargs prefetch --max-size 200GB
```

This will feed the accessions to the prefetch command, which will result in the raw sequences being downloaded from NCBI’s servers. You should see some progress messages as the files download.
 
Now that the SRA files should be downloaded to your scratch drive, we can extract the paired FASTQ files from the raw prefetched data. At this point it’s important to note what type of reads you are expecting. You’ll have to ensure that you get paired end reads from the SRA accessions where they are expected. The newer program called fasterq-dump appears to be aware of paired-end datasets, and splits them accordingly by default. 
Move into the directory containing all of the downloaded SRR files:
```
/scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/
```
And run the following command to extract the paired end FASTQ data from each accession. 

```
ls -1 | xargs -I'{}' fasterq-dump -t /scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/ -e 6 {}
```

By default, fasterq-dump uses 6 threads, but you can specify a different amount using the -e flag.

Using the -I'{}' flag of xargs, your accessions will be extracted serially, which is what we want in this case, in order to reduce the load on the server.
Note, you may need to increase the max limit size of the fastq files (20GB by default), depending on the size of your files. 

Now, you should be able to ls and see your brand new fastq files named with the SRRXXXXX_1.fastq for single end data, and SRRXXXXX_2.fastq for paired end data. 
```
ls -1 /scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/*fastq

SRR6656140.sra_1.fastq  
SRR6656155.sra_1.fastq  
SRR6656170.sra_1.fastq  
SRR6656185.sra_1.fastq  
SRR6656200.sra_1.fastq
```
#DONE

STEP-4: QC and mapping of SRAs

Now that the population level data has been downloaded and converted to paired-end fastq file format, the next step is to clean the reads of low quality base calls (<Q20) and remove any adapters present in the data (while these would likely be soft masked during the alignment step, best to just remove them here). This will all be accomplished using the trimgalor program, which will auto predict the adapter sequence based upon location and prevelance of nucleotides in the reads.

Copy / move the [trimgalor.sh](./trimgalore.sh) script to the following path, e.g.:
```
nano /scratch/bell/blackan/theta/GCF_002007445.1_ASM200744v2/sra/raw/trimgalore.sh
```
And then execute the SLURMM script as a job:



- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)

#### pipeline is tested until this step ####
QC of bam dataset
- qc_bams.sh & quantile_thresholds.R



Analyses
- angsd.sh



