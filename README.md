# Theta project

These scripts are developed for investigating genomic diversity and effective population size (θ = 4Neμ) in wildlife.

Jan 2021

The series of scripts will auto-create the following directory structure:
```

theta
  >species_x
      >species_x_ref
      >species_x_rm
      >species_x_gtf
      >sra
         >raw
         >cleaned
         >aligned
```
         
Once the folders are created, confirmed that all looks OK before proceeding through step 1-step X (NOTE, these scripts assume that a SLURMM scheduler is used; if a different scheduler is employed please modify script headers accordingly).

## Steps for analyzing each species

STEP 1: Download reference genome, repeater masker out file and annotations (if available):

a) For a current list of all mammalian refseq assemblies, go to: https://www.ncbi.nlm.nih.gov/datasets/genomes/?txid=40674&source_db=RefSeq
and locate your species of interest. 

b) Identify the best / most current assembly (based upon N50, assembly levlel, etc) and click on the refseq hyperlink beginning with "GCF_", which will link to the NCBI genome summary web-page. Once here, locate and navigate to the "FTP directory for RefSeq assembly" in the right hand column. This will pull up the File Transfer Protocal web page where the reference genome (e.g. GCF_002007445.1_ASM200744v2_genomic.fna.gz), repeat masker out (e.g. GCF_002007445.1_ASM200744v2_rm.out.gz) and annotations (e.g. GCF_002007445.1_ASM200744v2_genomic.gtf.gz) are available for download. For each file, copy the link address and paste each into a text editor in an open shell, such as:

```

#!/bin/sh -l
#SBATCH -A standby
#SBATCH --time=4:00:00
#SBATCH --job-name GCF_002007445.1_ASM200744v2_download

cd $SLURM_SUBMIT_DIR

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_genomic.fna.gz" \
-O GCF_002007445.1_ASM200744v2_ref/GCF_002007445.1_ASM200744v2_genomic.fna.gz -o log_ref 

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_rm.out.gz" \
-O GCF_002007445.1_ASM200744v2_rm/GCF_002007445.1_ASM200744v2_rm.out.gz -o log_rm

 wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/007/445/GCF_002007445.1_ASM200744v2/GCF_002007445.1_ASM200744v2_genomic.gtf.gz" \
 -O GCF_002007445.1_ASM200744v2_rm/GCF_002007445.1_ASM200744v2_genomic.gtf.gz -o log_gtf
 
 gzip -d */*gz

```

Save the file as "GCF_002007445.1_ASM200744v2_download.sh" within the "GCF_002007445.1_ASM200744v2" directory before running script

Once completed, check the log files to ensure that all fully downloaded (100%) and were decompressed!

```

cat log* | grep "100%" 
#There should be one line printed out for each file downloaded

ls */
#All files should be decompressed

```

STEP 2:

Now that the files have been download and decompressed, the reference needs to be QC'd using the following scripts:

  a) [qc_reference.sh](./qc_reference.sh) & [repeatmasker_names.R](./repeatmasker_names.R)
  
  b) [qc_reference_stats.sh](./qc_reference_stats.sh) & [qc_reference_stats.R](./qc_reference_stats.R)

Follow the pipeline below to qc reference assembly, qc read files, map reads, analyse data etc.

QC of reference assembly


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



