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
//currently it reads an input CSV, processes genome data, maps sequences, 
//then finally output genomic diversity metrics
========================================================================================
========================================================================================
*/

/*
========================================================================================
========================================================================================
//PROCESSES FOR THE PIPELINE
//PIPELINE EXECUTION CONTROL IN WORKFLOW AT THE END
========================================================================================
    SETUP A NEEDED DIRECTORY (e.g., /scratch/negishi/jeon96/theta/) 
    AND STORE NEEDED SCRIPTS IN IT
========================================================================================
*/

process step1{
    tag "$step1"
    clusterOptions '--job-name=Step1 -n 32 -N 1 -t 1-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL' 
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    cd /scratch/negishi/jeon96/theta/
    bash step1.sh ${genus_species} ${accession} ${pathway} ${assembly}
    """
}

process step2{
    tag "$step2"
    clusterOptions '--job-name=Step2 -n 32 -N 1 -t 1-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    bash step2.sh ${genus_species} ${accession} ${pathway} ${assembly}
    """
}

process step3{
    tag "$step3"
    clusterOptions '--job-name=Step3 -n 64 -N 1 -t 5-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly) 

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    bash step3.sh ${genus_species} ${accession} ${pathway} ${assembly}
    """
}

process step4{
    tag "$step4"
    clusterOptions '--job-name=Step4 -n 64 -N 1 -t 1-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly) 

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    bash step4.sh ${genus_species} ${accession} ${pathway} ${assembly}
    """
}

process step5{
    tag "$step5"
    clusterOptions '--job-name=Step5 -n 64 -N 1 -t 5-00:00:00 -A fnrdewoody'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(genus_species), val(accession), val(pathway), val(assembly) 

    output:
    tuple val(genus_species), val(accession), val(pathway), val(assembly)

    script:
    """
    bash step5.sh ${genus_species} ${accession} ${pathway} ${assembly}
    """
}

/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    //input channel for a CSV that could contain parameters for Theta workflow
    print "Welcome to this CSV workflow to estimate several GD metrics"
    print 'You will be working on the following files: '+params.csv
    Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> [row.genus_species, row.accession, row.pathway, row.assembly] } \
        | view() \
        //| set {csv} \
	| step1 \
	| step2 \
	| step3 \
        | step4 \
        | step5
}
