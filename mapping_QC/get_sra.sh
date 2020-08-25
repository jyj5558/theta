#!/bin/bash
#SBATCH --job-name=sra
#SBATCH -A fnrgenetics
#SBATCH -t 300:00:00 
#SBATCH -N 1 
#SBATCH -n 1

module load bioinfo

# cd $SLURM_SUBMIT_DIR

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1013736/ERR1013736_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1013736/ERR1013736_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/007/ERR1013737/ERR1013737_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/007/ERR1013737/ERR1013737_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1013738/ERR1013738_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1013738/ERR1013738_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013739/ERR1013739_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013739/ERR1013739_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR1013740/ERR1013740_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR1013740/ERR1013740_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/001/ERR1013741/ERR1013741_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/001/ERR1013741/ERR1013741_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1013742/ERR1013742_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/002/ERR1013742/ERR1013742_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1013743/ERR1013743_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/003/ERR1013743/ERR1013743_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1013744/ERR1013744_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/004/ERR1013744/ERR1013744_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/005/ERR1013745/ERR1013745_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/005/ERR1013745/ERR1013745_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1013746/ERR1013746_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/006/ERR1013746/ERR1013746_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/007/ERR1013747/ERR1013747_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/007/ERR1013747/ERR1013747_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1013748/ERR1013748_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/008/ERR1013748/ERR1013748_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013749/ERR1013749_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013749/ERR1013749_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR1013750/ERR1013750_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/000/ERR1013750/ERR1013750_1.fastq.gz
