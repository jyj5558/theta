# Theta project
Jan 2021

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife populations. NOTE, these scripts assume that a SLURM scheduler is used; if a different scheduler is employed (e.g., SGE) please modify script headers and module loads accordingly).

## Steps for analyzing each species (one at a time)

STEP 1: Start an interactive session:

Begin an interactive session to move to a back-end machine:
```
#Start an interactive session on standby
sinteractive -t 4:00:00 --mem=10G
```

STEP 2: Create folder structure and download reference genome, repeater masker out file and annotations (if available). NOTE, below the Cheetah genome is used as an example only. User's will need to designate the genome assembly for the species they are targetting. 

a) For a current list of all mammalian refseq assemblies, go to: https://www.ncbi.nlm.nih.gov/datasets/genomes/?txid=40674&source_db=RefSeq
and locate your species of interest. 

b) Identify the best / most current assembly (based upon N50, assembly levlel, etc) and click on the refseq hyperlink beginning with "GCF_", which will link to the NCBI genome summary web-page. 

Once here, locate and navigate to the "FTP directory for RefSeq assembly" in the right hand column. This will pull up the File Transfer Protocal web page where the reference genome (e.g. GCF_002007445.1_ASM200744v2_genomic.fna.gz), repeat masker out (e.g. GCF_002007445.1_ASM200744v2_rm.out.gz) and annotations (e.g. GCF_002007445.1_ASM200744v2_genomic.gtf.gz) are available for download. 

Copy the name of the reference genome build, and the link addresses to all three files (if available)

Now that you have this information recorded, create the [setup.sh](./setup.sh) script in the ./theta folder and make the script executable:
```
chmod +x ./setup.sh
```

Now, run the [setup.sh](./setup.sh) script in the ./theta folder to 1) auto-create the proper directory structure and b) download and decompress all files. 

The user will be prompted to enter the following information, which will create the necessary folder structure, download all files and decompress them. If the user species does not have one of these listed files, just enter NA when prompted).

For example:

```
./setup.sh

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

Once completed, check the standard out to ensure that all fully downloaded (100%), decompressed, named correctly and are not empty. For example, the user should see something like this printed out on the screen when running the above script:

```
Please confirm that all downloaded files were fully downloaded 100%
 15450K .......... .                                          100%  375M=0.3s
733150K .......... .......... .......... ....                 100%  278M=11s
128950K .....                                                 100%  379M=1.9s
```
And confirm that none of the files are empty (based upon file size) and named accordingly:

```
Are all downloaded files named correctly and are not empty?
GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_gtf:
total 400M
-rw-r--r-- 1 blackan student 500M Oct 29 03:52 GCF_009873245.2_mBalMus1.pri.v3.gtf

GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_ref:
total 2.2G
-rw-r--r-- 1 blackan student 2.3G Oct 29 03:52 GCF_009873245.2_mBalMus1.pri.v3.genomic.fna

GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_rm:
total 464M
-rw-r--r-- 1 blackan student 544M Oct 29 03:52 GCF_009873245.2_mBalMus1.pri.v3.rm.out
```
Finally, the header of the mitochondrial sequence (if present) will be identified and used to remove this sequence from the reference genome. The script will print out information on this step, so please check for any memory issues or error messages:

```
Now, the reference will be searched for a mitochondrion sequence
>NC_001601.1 Balaenoptera musculus mitochondrion, complete genome
 
If a sequence was printed out, this will now be removed from the reference genome
java -ea -Xmx8191m -cp /group/bioinfo/apps/apps/BBMap-37.93/current/ driver.FilterReadsByName include=f in=./GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_ref/GCF_009873245.2_mBalMus1.pri.v3.genomic.fna out=./GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_ref/tmp.fasta names=mitochondrion ow=t substring=t
Executing driver.FilterReadsByName [include=f, in=./GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_ref/GCF_009873245.2_mBalMus1.pri.v3.genomic.fna, out=./GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_ref/tmp.fasta, names=mitochondrion, ow=t, substring=t]

Input is being processed as unpaired
Time:               9.713 seconds.
Reads Processed:    106 	0.01k reads/sec
Bases Processed:    2374868943 	244.50m bases/sec
Reads Out:          105
Bases Out:          2374852541
```

STEP 2 (*UNTESTED*):

Now that the files have been download and decompressed, the reference needs to be QC'd using the following scripts:

  a) [qc_reference.sh](./qc_reference.sh) & [repeatmasker_names.R](./repeatmasker_names.R)
  
  b) [qc_reference_stats.sh](./qc_reference_stats.sh) & [qc_reference_stats.R](./qc_reference_stats.R)


STEP 3: Locate and download SRA files for population level reads:

Locate SRA files for a given species:
NEEDS TO BE WRITTEN!

Download SRA files. 
Note, as per opinion piece linked [here](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13309?campaign=wolearlyview) SRA files need to be categorized according to the sequencing platform they were generated on (i.e., 2- or 4-channel chemistry) to avoid biases in generating theta estimates. 

Therefore, please create two lists (if applicable): 1) a comma seperated list of SRA files generated using 2-channel chemistry and 2) a comma seperated list of SRA files generated using 4-channel chemistry. 

For example, 

```
#Two-channel chemistry SRA names
ERR1843087,ERR1843088,ERR1843089
```

```
#Four-channel chemistry SRA names
SRR10066840,SRR10066841,SRR10066842
```
Now, the user will need to download the SRA files using an interactive session:
```
sinteractive -A standby -t 4:00:00 -n 6 --mem=10G
```

Within the main directory (e.g., theta)
```
cd /scratch/bell/dewoody/theta
```

The user can run the [download_sra_2.sh](./download_sra_2.sh)  script, such as:
```
chmod +x ./download_sra_2.sh
./download_sra_2.sh
```
The user will be prompted to paste in the name of the genome assembly, such as:

```
Hello, what is the name of the accession for this species genome?
GCF_009873245.2_mBalMus1.pri.v3
```
This will change the working directory to this species download folder.
The user will then be prompted to enter the CSV list of SRA accessions:
```
Please specify the name of the file containing sra names, this should be located in your pwd directory for GCF_009873245.2_mBalMus1.pri.v3
ERR1843087,ERR1843088,ERR1843089
```
This script will then download all SRA files and dump them into paired-end fastq file format for TWO-chemistry platforms.

Repeat steps above if the user also has 4-channel SRA files to download, with the [download_sra_4.sh](./download_sra_4.sh) 

Now, there should be fully downloaded fastq files named with the [S,E]RRXXXXX_1.fastq for single end data, and [S,E]RRXXXXX_2.fastq for paired end data within the ./two_ch/raw and/or ./four_ch/raw directories. 

Note, make sure that there are no unresolved downloads during either the SRA download or fastq!


STEP-4: QC and mapping of SRAs

Now that the population level data has been downloaded and converted to paired-end fastq file format, the next step is to clean the reads of low quality base calls (<Q20) and remove any adapters present in the data (while these would likely be soft masked during the alignment step, best to just remove them here). This will all be accomplished using the trimgalor program, which will auto predict the adapter sequence based upon location and prevelance of nucleotides in the reads.

Copy / move the [cleaning.sh](./cleaning.sh) script to the following path, e.g.:
```
/scratch/bell/dewoody/theta/GCF_009873245.2_mBalMus1.pri.v3/sra/raw
```
And then execute the SLURMM script as a job:

```
sbatch cleaning.sh
```

Once the job is finished, check the output of the slurm job to get a snapshot of the cleaning metrics. 
User should be in the 'raw' directory and run the following command:
```
ls -1 *fastq > tmp.1 ; grep "Total written (filtered):" slurm* > tmp.2 ; paste tmp.1 tmp.2 ; rm tmp.*
```
Which should print out someting like this:

```
ERR1843087.sra_1.fastq	Total written (filtered):     92,177,528 bp (95.8%)
ERR1843087.sra_2.fastq	Total written (filtered):     83,786,746 bp (87.5%)
. . .
```
#DONE

#### pipeline is tested until this step ####




- mapping.sh
- realignment_single.sh (for individuals with a single SRA)
- realignment_multiple.sh (for individuals with multiple SRAs)


QC of bam dataset
- qc_bams.sh & quantile_thresholds.R



Analyses
- angsd.sh



