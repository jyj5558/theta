#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//where are the inputs located
//path to CSV file
params.csv = "$CLUSTER_SCRATCH"

//option to change from default directory to path of choice
params.savepath = "$CLUSTER_SCRATCH"

//user email information
params.email = "" //your email address

/*
========================================================================================
========================================================================================
//A TEMPLATE FOR A NEXTFLOW WORKFLOW FOR USE ON PURDUE RCAC CLUSTER
//DESCRIPTION: WHAT DOES THIS WORKFLOW DO?
//currently it just reads an input CSV...stay tuned for more
========================================================================================
========================================================================================
*/

/*
========================================================================================
========================================================================================
//PROCESSES FOR THE PIPELINE
//PIPELINE EXECUTION CONTROL IN WORKFLOW AT THE END
========================================================================================
    SETUP A NEEDED DIRECTORY
========================================================================================
*/

process step1{
    tag "$step1"
    clusterOptions '--job-name=Step1 -n 32 -t 1-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL' 
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    module load biocontainers
    module load bioawk
    module load samtools
    module load bbmap

    ####notes and usage#####
    #This script downloads reference, and repeat 
    #If a masked genome isn't available (i.e. rm.out), script will create one using the mammal database
    #Script was written to be run with SLURM job scheduler
    #
    #User will need to input (paste) information for the following variables below:
    #
    #genus-species: this is used in the directory naming, to make browsing a bit more doable for us humans
    #accession: this is also used in the directory, to keep multiple reference assemblies separate
    #pathway: include full NCBI url to FTP site (containing directory)                                          
    #assembly: name of assembly
    #n: the number of cpus you allocated in the SBATCH command above
    #
    #Example of defined variables below 
    #genus_species=Balaenoptera-musculus
    #accession=GCF_009873245.2
    #pathway=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/
    #assembly=mBalMus1.pri.v3
    #n=32
    ####end of notes and usage####


    ########################
    #DO NOT EDIT BELOW CODE#
    ########################

    #Move to the scratch space
    cd $CLUSTER_SCRATCH

    ####create directories and download reference genome, repeat masker, and annotation####
    mkdir -p ./$genus_species/${accession}_ref
    mkdir ./$genus_species/${accession}_rm
    mkdir ./$genus_species/${accession}_gtf
    cd $genus_species

    #Download reference genome
    cd ${accession}_ref
    wget ${pathway}${accession}_${assembly}_genomic.fna.gz 
    gunzip ${accession}_${assembly}_genomic.fna.gz
    cp ${accession}_${assembly}_genomic.fna original.fa # keep a copy of the original reference
    cd ../

    #Download repeatmasker file (if available)
    cd ${accession}_rm
    wget ${pathway}${accession}_${assembly}_rm.out.gz 
    gunzip ${accession}_${assembly}_rm.out.gz
    cp ${accession}_${assembly}_rm.out rm.out # keep a copy of the original repeatmasker
    cd ../

    #Print out file sizes for checking later
    ls -lh ${accession}* > download_log

    #Search for and remove mito seq in ref genome. Will only work if marked in assembly!
    grep "mitochondrion" ${accession}_ref/original.fa | cut -f 1 > mito_header.txt #If no mitochondrial sequence present in reference, will be blank. Otherwise contain header
    filterbyname.sh include=f in=${accession}_ref/original.fa out=${accession}_ref/original.tmp.fa names="mitochondrion" ow=t substring=t
    rm ${accession}_ref/original.fa
    mv ${accession}_ref/original.tmp.fa ${accession}_ref/original.fa      

    ###prep reference genome for mapping####
    #Reduce fasta header length
    reformat.sh in=${accession}_ref/original.fa out=${accession}_ref/new.fa trd=t -Xmx20g overwrite=T
    #Sort by length
    sortbyname.sh in=${accession}_ref/new.fa out=${accession}_ref/ref.fa -Xmx20g length descending overwrite=T
    #Remove sequences smaller that 100kb prior to any repeatmasking
    bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' ${accession}_ref/ref.fa > ${accession}_ref/ref_100kb.fa
    rm ${accession}_ref/new.fa

    #Index ref.fa and ref_100kb.fa for step3, step4, and step5 
    samtools faidx ${accession}_ref/ref_100kb.fa
    samtools faidx ${accession}_ref/ref.fa

    #Prep repeatmasked file for later processing, create a rm.out if one is not available. 
    #Move into rm directory
    cd ${accession}_rm/ 

    FILE1=$"rm.out"
    if [ -s $FILE1 ]
    then
	    # parse rm.out out to create bed coordinates
	    cat rm.out |tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' > repeats.bed # make bed file
    else	
	    #If no rm.out file is available, run RepeatMasker. Note, samtools conflict so had to purge first
	    module --force purge
	    module load biocontainers/default
	    module load repeatmasker
	    RepeatMasker -pa $n -a -qq -species mammals -dir . ../${accession}_ref/ref_100kb.fa 
	    cat ref_100kb.fa.out  | tail -n +4 | awk '{print $5,$6,$7,$11}' | sed 's/ /\t/g' > repeats.bed 
    fi

    #Move back to species/accession directory
    cd ../ 

    #END
    """
}

process step2{
    tag "$step2"
    clusterOptions '--job-name=Step2 -n 32 -t 5-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    module load biocontainers
    module load trim-galore
    module load cutadapt
    module load fastqc
    module load sra-tools

    ####notes and usage####
    #
    #Edit the step2_SRA_download_clean.sh script, to add target species "genus-species" and the number of cpus "n" you allocated in the SBATCH command above
    #script was written to be run with SLURM job scheduler
    #
    #Example:
    #genus_species=Marmota-marmota-marmota
    #n=32
    #theta git repository should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/abc/xyz/)
    #Genus-species: this is used in the directory naming to make browsing a bit more doable for us humans
    #
    ####end of notes and usage####

    ########################
    #DO NOT EDIT BELOW CODE#
    ########################

    #Extract relevant info from metadata file in theta directory for target species and save in species folder:

    cd $CLUSTER_SCRATCH/theta/${genus_species}/
    cat $CLUSTER_SCRATCH/theta/SRA_metadata/${genus_species}.txt | sed 's/ /_/g' > ${genus_species}_SRA.txt


    ####create directories and download SRAs####
    mkdir -p ./sra/raw
    mkdir ./sra/cleaned
    mkdir ./sra/aligned

    cd ./sra/raw/

    cat ../../${genus_species}_SRA.txt | cut -f 1 | tail -n +2 | while read g
    do
    mkdir ${g}
    cd ${g}
    prefetch --max-size 500GB -O ./ ${g}
    fasterq-dump -e $n --progress ${g}.sra
    find . -name '*.fastq' -exec mv {} ../ \;
    cd ../
    rm -r ${g}

    fastqc ${g}_1.fastq --extract --quiet
    fastqc ${g}_2.fastq --extract --quiet
    rm ${g}_1_fastqc.zip
    rm ${g}_2_fastqc.zip

    #Merge output from fastqc and check for FAILs
    cat ${g}_1_fastqc/summary.txt ${g}_2_fastqc/summary.txt > ${g}_fastqc_summary.txt
    FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
    echo "raw"
    echo "$FAIL"
    rm -r ${g}_?_fastqc* 

    #Trim adapters
    trim_galore --stringency 1 --length 30 --quality 20 --fastqc_args "--nogroup" -o ../cleaned --paired ${g}_1.fastq ${g}_2.fastq

    #Check quality of trimmed reads
    cd ../cleaned
    cat ${g}_1.fastq_trimming_report.txt ${g}_2.fastq_trimming_report.txt > ${g}_fastqc_summary.txt
    rm ${g}_1_val_1_fastqc.zip
    rm ${g}_2_val_2_fastqc.zip
    FAIL=$(grep "FAIL" ${g}_fastqc_summary.txt)
    echo "cleaned"
    echo "$FAIL"
    rm ${g}_fastqc_summary.txt
    rm ${g}_1_val_1_fastqc.html
    rm ${g}_2_val_2_fastqc.html
    cd ../raw
    done
    #END
    """
}

process step3{
    tag "$step3"
    clusterOptions '--job-name=Step3 -n 64 -t 14-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly) 

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    module load biocontainers
    module load bwa
    module load samtools
    module load picard
    module load gatk
    module load bamtools

    ####notes and usage####
    #
    # this script maps read files to the QC reference assembly					#
    # the following files to use in downstream analyses are created:	  	#
    # *.sam file for each SRR*												#
    # script was written to be run with SLURM job scheduler
    #
    # add target species "genus-species"
    #Example:
    #genus_species=Marmota-marmota-marmota
    #theta git should be cloned at $CLUSTER_SCRATCH (e.g., /scratch/abc/xyz/)
    ####end of notes and usage####

    ########################
    #DO NOT EDIT BELOW CODE#
    ########################

    #Make a new directory to house processed alignment files
    mkdir $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/

    #Move to reference

    #Index/create dict files
    cd $CLUSTER_SCRATCH/theta/${genus_species}/${accession}_ref/

    #Create dictionary for realignment
    PicardCommandLine CreateSequenceDictionary reference=ref.fa output=ref.dict

    #Move to cleaned fastq files in preparation for alignment
    cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/cleaned/

    #Index reference 
    bwa index -a bwtsw ../../${accession}_ref/ref.fa

    #Capture all cleaned fastq files with variable
    ls -1 *.fq | sed "s/_[1-2]_val_[1-2].fq//g" | uniq > cleaned_sralist
    for i in `cat cleaned_sralist`
    do

    #Perform alignment using whole CPUs and bwa mem algorithm
    bwa mem -t 64 -M -R "@RG\tID:group1\tSM:${i}\tPL:illumina\tLB:lib1\tPU:unit1" ../../${accession}_ref/ref.fa  ${i}_1_val_1.fq ${i}_2_val_2.fq > ../aligned/${i}.sam

    #Move to the directory containing the alignment files
    cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/aligned/

    #Validated *sam files should produce a summary output file that contains no errors
    PicardCommandLine ValidateSamFile I=${i}.sam MODE=SUMMARY O=${i}.sam.txt

    #Check for errors in sam file
    noerror=$(grep "No errors found" ${i}.sam.txt)
    [[ ! "$noerror" ]] && echo "$i samfile not OK" || echo "$i samfile OK"

    #Sort sam file based upon read coordinate
    mkdir -p ./tmp/
    PicardCommandLine SortSam TMP_DIR=$CLUSTER_SCRATCH/theta/${genus_species}/sra/aligned/tmp INPUT=${i}.sam OUTPUT=${i}_sorted.bam SORT_ORDER=coordinate

    #Mark PCR duplicates without removing them
    PicardCommandLine MarkDuplicates INPUT=${i}_sorted.bam OUTPUT=./${i}_marked.bam METRICS_FILE=${i}_metrics.txt
    PicardCommandLine BuildBamIndex INPUT=./${i}_marked.bam

    # local realignment of reads
    GenomeAnalysisTK -nt 64 -T RealignerTargetCreator -R ../../${accession}_ref/ref.fa -I ${i}_marked.bam -o ${i}_forIndelRealigner.intervals

    #Realign with established intervals
    GenomeAnalysisTK -T IndelRealigner -R ../../${accession}_ref/ref.fa -I ${i}_marked.bam -targetIntervals ${i}_forIndelRealigner.intervals -o ../final_bams/${i}.bam

    cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/

    #Get some summar stats on bam files
    echo "$i mapping rate is" > ./${i}_mapping.txt
    samtools flagstat ./${i}.bam >> ./${i}_mapping.txt
    echo "$i depth is" > ./${i}_depth.txt
    samtools depth -a ./${i}.bam | awk '{c++;s+=$3}END{print s/c}' >> ./${i}_depth.txt
    echo "$i breadth is" > ./${i}_breadth.txt
    samtools depth -a ./${i}.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> ./${i}_breadth.txt

    cd $CLUSTER_SCRATCH/theta/${genus_species}/sra/cleaned/
    done

    # END
    """
}

process step4{
    tag "$step4"
    clusterOptions '--job-name=Step4 -n 64 -t 14-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """  
    module load biocontainers
    module load bioawk
    module load cmake/3.9.4
    module load bedtools
    module load bbmap
    module load R
    module load bedops
    export PATH=$PATH:~/genmap-build/bin

    ####notes####
    #
    # This script estimates mappability and finds all of the short scaffolds and then identifyies repeats. 
    # the package GenMap should be downloaded in advance at the home directory
    #
    # The output files include: 								
    # ok.bed (regions to analyze in angsd etc)		
    # map_repeat_summary.txt (summary of ref quality)							
    # script was written to be run with SLURM job scheduler
    #
    ####usage####
    #User will need to input (paste) information for the following variables below:
    #
    #genus_species: this is used in the directory naming to make browsing a bit more doable for us humans
    #accession: this is also used in the directory, to keep multiple reference assemblies separate
    #n: the number of cpus you allocated in the SBATCH command above
    #
    #Example of defined variables below 
    #
    #genus_species=Balaenoptera-musculus
    #accession=GCF_009873245.2
    #n=32
    ####end of notes and usage####

    ########################
    #DO NOT EDIT BELOW CODE#
    ########################

    #Move to species/accession directory

    cd $CLUSTER_SCRATCH/theta/${genus_species}/

    ####Assess mappability of reference####
    genmap index -F ${accession}_ref/ref_100kb.fa -I index -S 50 # build an index 

    #Compute mappability, k = kmer of 100bp, E = # two mismatches
    mkdir mappability
    genmap map -K 100 -E 2 -T $n -I index -O mappability -t -w -bg                

    #Sort bed 
    sortBed -i ${accession}_rm/repeats.bed > ${accession}_rm/repeats_sorted.bed 

    #Make ref.genome
    awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ${accession}_ref/ref_100kb.fa.fai > ${accession}_ref/ref.genome 

    #Sort genome file
    awk '{print $1, $2, $2}' ${accession}_ref/ref.genome > ${accession}_ref/ref2.genome
    sed -i 's/ /\t/g' ${accession}_ref/ref2.genome
    sortBed -i ${accession}_ref/ref2.genome > ${accession}_ref/ref3.genome
    awk '{print $1, $2 }' ${accession}_ref/ref3.genome > ${accession}_ref/ref_sorted.genome
    sed -i 's/ /\t/g' ${accession}_ref/ref_sorted.genome
    rm ${accession}_ref/ref.genome
    rm ${accession}_ref/ref2.genome
    rm ${accession}_ref/ref3.genome

    #Find nonrepeat regions
    bedtools complement -i ${accession}_rm/repeats_sorted.bed -g ${accession}_ref/ref_sorted.genome > ${accession}_rm/nonrepeat.bed

    #Clean mappability file, remove sites with <1 mappability                                                    
    awk '$4 == 1' mappability/ref_100kb.genmap.bedgraph > mappability/map.bed                                           
    awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' mappability/map.bed > mappability/mappability.bed

    #Sort mappability 
    sortBed -i mappability/mappability.bed > mappability/mappability2.bed
    sed -i 's/ /\t/g' mappability/mappability2.bed

    #Only include sites that are nonrepeats and have mappability ==1
    bedtools subtract -a mappability/mappability2.bed -b ${accession}_rm/repeats_sorted.bed > mappability/map_nonreapeat.bed

    #Sort file -- by chr then site
    bedtools sort -i mappability/map_nonreapeat.bed > mappability/filter_sorted.bed

    #Merge overlapping regions
    bedtools merge -i mappability/filter_sorted.bed > mappability/merged.bed

    #Make bed file with the 100k and merged.bed (no repeats, mappability =1 sites) from below
    awk '{ print $1, $2, $2 }' ./${accession}_ref/ref_100kb.fa.fai > ./${accession}_ref/ref_100kb.info

    #Replace column 2 with zeros
    awk '$2="0"' ./${accession}_ref/ref_100kb.info > ./${accession}_ref/ref_100kb.bed

    #Make tab delimited
    sed -i 's/ /\t/g' ./${accession}_ref/ref_100kb.bed

    #Only include scaffolds in merged.bed if they are in ref_100kb.bed
    bedtools intersect -a ./${accession}_ref/ref_100kb.bed -b ./mappability/merged.bed > ok.bed	

    #Make chrs.txt
    cut -f 1 ${accession}_ref/ref_100kb.bed > ${accession}_ref/chrs.txt
	
    #Remove excess files
    rm mappability/map.bed
    rm mappability/filter_sorted.bed
    rm mappability/mappability2.bed
    rm ${accession}_ref/ref_sorted.genome

    #Output some QC stats
    cd $CLUSTER_SCRATCH/theta/source/
    Rscript qc_reference_stats.R --args $CLUSTER_SCRATCH/theta/${genus_species}/ ${genus_species} ${accession} 
    cd $CLUSTER_SCRATCH/theta/${genus_species}/

    map=$(sed -n '1p' okmap.txt)
    norepeats=$(sed -n '1p' norepeat.txt)
    okbed=$(sed -n '1p' okbed.txt)

    echo -e "${genus_species}\t ${accession}\t ${map}\t ${norepeats}\t ${okbed}" >> map_repeat_summary.txt

    #Make list of the bamfiles and index each file in theta directory for step5 in advance
    # if you remove some individuals depending on their mapping rate, depth, or breadth, edit this file accordingly before step5
    mkdir $CLUSTER_SCRATCH/theta/${genus_species}/theta/
    cd $CLUSTER_SCRATCH/theta/${genus_species}/theta/
    ls $CLUSTER_SCRATCH/theta/${genus_species}/sra/final_bams/*.bam > ./bam.filelist

    # END
    """
}

process step5{
    tag "$step5"
    clusterOptions '--job-name=Step5 -n 64 -t 14-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly) 

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """  
    module --force purge
    module load biocontainers
    module load angsd
    module load r
    module load bcftools
    module load htslib

    ####notes and usage####
    #
    # this script identifies / filteres for usable sites and estimates theta#
    # the following result files are created:	
    # Thetas_${genus_species}.txt, IndHet_${genus_species}.txt, PopHet_${genus_species}.txt,
    # ROH_${genus_species}.txt
    # script was written to be run with SLURM job scheduler
    #
    ####end of notes and usage####

    ########################
    #DO NOT EDIT BELOW CODE#
    ########################

    #Absolute paths of parent and theta directory of genus-species shorten commands
    PD=$CLUSTER_SCRATCH/theta/${genus_species}
    THETA=$CLUSTER_SCRATCH/theta/${genus_species}/theta

    #Move to directory which will house angsd/theta files
    mkdir -p $THETA
    cd $THETA

    #Designate min number of individuals, set to total numb of bam-filed individuals divided by two
    MIND=$((`wc -l < ./bam.filelist` / 2))

    #Convert bed file to angsd format
    awk '{print $1"\t"$2+1"\t"$3}' $PD/ok.bed > ./angsd.file

    #Index file
    angsd sites index ./angsd.file

    #Estimate GL
    echo "Genotype likelihood estimation started"
    angsd -bam ./bam.filelist -ref $PD/${accession}_ref/ref.fa -anc $PD/${accession}_ref/ref.fa -rf $PD/${accession}_ref/chrs.txt -sites ./angsd.file \
    -dosaf 1 -GL 2 -doMajorMinor 1 \
    -minInd $MIND -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 \
    -out out -P 64 
    echo "Genotype likelihood estimation done"
    #Obtain ML estimate of SFS using the folded realSFS
    echo "SFS estimation started"
    realSFS -P 64 out.saf.idx  -fold 1 > out.sfs
    echo "SFS estimation done"

    #Calculate theta for each site
    echo "Theta estimation started"
    realSFS saf2theta out.saf.idx -sfs out.sfs -outname out -P 64 

    #Estimate 
    thetaStat print out.thetas.idx > out.thetas_persite.txt
    echo "per-site Theta file created"

    #Sliding window estimate
    thetaStat do_stat out.thetas.idx -win 50000 -step 50000 -outnames theta.thetasWin.gz
    echo "window-based Theta file created"

    #Column 4 has Wattersons, column 5 has Nucleotide diversity, column 9 has Tajima's D, and column 10 has Fu & Li's Fs
    awk '{print $4,$5,$9,$10,$14}' theta.thetasWin.gz.pestPG > Thetas

    #Get mean
    sumSites=$(awk 'BEGIN{s=0;}{s=s+$5;}END{print s;}' Thetas)
    sumW=$(awk 'BEGIN{w=0;}{w=w+$1;}END{print w;}' Thetas)
    sumN=$(awk 'BEGIN{n=0;}{n=n+$2;}END{print n;}' Thetas)
    meanW=$(awk "BEGIN {print $sumW/$sumSites}")
    meanN=$(awk "BEGIN {print $sumN/$sumSites}")
    meanD=$(awk 'BEGIN{d=0;}{d=d+$3;}END{print d/NR;}' Thetas)
    meanF=$(awk 'BEGIN{f=0;}{f=f+$4;}END{print f/NR;}' Thetas)

    echo -e "Genus_Species\tNucleotide_Diversity\tWatterson_Theta\tTajima_D\tFu&Li_F\n${genus_species}\t$meanN\t$meanW\t$meanD\t$meanF\t" >> Thetas_${genus_species}.txt
    echo "Thetas (Pi, Theta, D, F) .txt created"

    #Heterozygosity for each individual

    mkdir ./HET
    OUTDIR='HET'

    cat ./bam.filelist | sed 's/\//\'$'\t/g' | cut -f 9 | sed 's/.bam//g' | while read -r LINE

    do
    echo "${LINE} heterozygosity estimation started"

    angsd -i ../sra/final_bams/${LINE}.bam -ref $PD/${accession}_ref/ref.fa -anc $PD/${accession}_ref/ref.fa -rf $PD/${accession}_ref/chrs.txt -sites ./angsd.file \
    -dosaf 1 -GL 2 -doMajorMinor 1 \
    -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1 -baq 2 -doCounts 1 -setMinDepthInd 5 \
    -out ${OUTDIR}/${LINE} -P 64
 
    realSFS -P 64 -fold 1 ${OUTDIR}/${LINE}.saf.idx > ${OUTDIR}/${LINE}_est.ml

    cd ${OUTDIR}
    Rscript -e 'args<-commandArgs(TRUE); LINE<-args[1]; a<-scan(paste(LINE,"est.ml", sep="_")); a[2]/sum(a)' ${LINE} >>  ../Het
    echo "${LINE} heterozygosity estimation done"
    cd ../
    done

    cat ./Het | cut -d " " -f 2 | tr '\n' ',' > ./IndHet_${genus_species}.txt

    #Heterozygosity for population
    #Get mean
    meanH=$(awk 'BEGIN{s=0;}{s=s+$2;}END{print s/NR;}' Het)

    #Get SD
    sdH=$(awk '{delta = $2 - avg; avg += delta / NR; \
    meanH2 += delta * ($2 - avg); } END { print sqrt(meanH2 / NR); }' Het)

    #Print to file
    echo -e "$PWD\t $meanH\t $sdH" \
    >> PopHet_${genus_species}.txt
    echo "Population heterozygosity .txt created"

    #ROHs

    echo "File conversion to bcf started"
    angsd -b ./bam.filelist -ref $PD/${accession}_ref/ref.fa -rf $PD/${accession}_ref/chrs.txt -sites ./angsd.file \
    -dobcf 1 -gl 2 -dopost 1 -domajorminor 1 -domaf 1 \
    -minMapQ 30 -minQ 30 -only_proper_pairs 1 -remove_bads 1 -uniqueOnly 1  -baq 2 -snp_pval 1e-6 -P 64
    echo "bcf file created"

    echo "ROH estimation started"
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' angsdput.bcf | bgzip -c > ${genus_species}.freqs.tab.gz
    tabix -s1 -b2 -e2 ${genus_species}.freqs.tab.gz
    bcftools roh --AF-file ${genus_species}.freqs.tab.gz --output ROH_${genus_species}_PLraw.txt --threads 64 angsdput.bcf
    echo "ROH raw files created"

    echo "ROH raw file parsing started"
    python3 $CLUSTER_SCRATCH/theta/ROHparser.py ${genus_species} ${accession}
    # END
    """
}

/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    //input channel for a CSV that could contain parameters for Jong Yoon's workflow
    print "Welcome to this CSV workflow to estimate several GD metrics"
    print 'You will be working on the following files: '+params.csv
    Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> [row.genus_species, row.accession, row.pathway, row.assembly] } \
        | view() \
        | set {csv}
    step1(csv) | set{step1Out}
    step2(csv) | set{step2Out} 
    step3In = step1Out.mix(step2Out) | uniq() | view() 
    step3(step3In) | set{step3Out} # use view; # test each piece first 
    step4(step1Out) | set{step4Out}
    step5In = step3Out.mix(step4Out) | uniq() | view()
    step5(csv)
}
