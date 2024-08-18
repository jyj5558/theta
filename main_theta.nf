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
    clusterOptions '--job-name=Step1 -n 32 -N 1 -t 2-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL' 
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
    clusterOptions '--job-name=Step2 -n 32 -N 1 -t 5-00:00:00 -A fnrdewoody'
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
    clusterOptions '--job-name=Step3 -n 64 -N 1 -t 14-00:00:00 -A fnrdewoody'
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
    clusterOptions '--job-name=Step4 -n 64 -N 1 -t 4:00:00 -A standby'
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
    clusterOptions '--job-name=Step5 -n 64 -N 1 -t 14-00:00:00 -A fnrdewoody'
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
        | set {csv}
    step1(csv) | set {step1_out}
    step2(csv) | set {step2_out} //run step1 and step2 in parallel    
    step1_out.mix(step2_out) | collect(flat: false) | flatMap | unique | set {step3_in} //check and collect step1 and step2 outputs
    step3_in.view()
    step3(step3_in) | set {step3_out}
    step4(step1_out) | set {step4_out} //step4 just processes step1 outputs
    step3_out.mix(step4_out) | collect(flat: false) | flatMap | unique | set {step5_in} //check and collect step3 and step4 outputs
    step5_in.view()
    step5(step5_in)
    
    //if above does not work, use below pipeline (without double slashes) after "view() \" above.
    
	//| step1 \
	//| step2 \
	//| step3 \
        //| step4 \
        //| step5    
    
}

