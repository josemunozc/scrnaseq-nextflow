#!/usr/bin/env nextflow

 println """\
===================================
 S C R N A S E Q - P I P E L I N E
===================================
General Information:
--------------------
  Profile(s):                 ${workflow.profile}
Input Parameters:
-----------------
  Base Directory:             ${params.projectDir}
  FastQC path:                ${params.dataDir}

"""

Channel.fromPath("${params.dataDir}/samplesheet_1kmouse.csv")
       .splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .into { sampleNames }

// copy genome files

process copy_genome {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'htc'
    time params.copyGenomeJobLength

    tag "Copying genome file into resources/"
    publishDir path:{params.resourcesDir},mode: 'copy'

    output:
    file("${params.genomeName}") into genome_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    cp ${params.genome} .

    sleep ${params.sleepTimeEnd}
    """
}

process copy_gtf {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'htc'
    time params.copyGTFJobLength

    tag "Copying gtf file into resources/"
    publishDir path:{params.resourcesDir},mode: 'copy'

    output:
    file("*") into gtf_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    cp ${params.gtf} .

    sleep ${params.sleepTimeEnd}
    """
}

process make_reference {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1  
    executor 'slurm'
    queue 'htc'
    time params.mkrefJobLength
    memory params.mkrefMemory

    tag "Making reference files"
    publishDir path:{params.resourcesDir},mode: 'copy'

    input:
    file genome from genome_ch
    file gtf from gtf_ch

    output:
    file("${params.organism}") into mkref_ch

    script:
    """

    ${cellrangerPath}/cellranger mkref --genome=${params.organism} --fasta=${genome} --genes=${gtf}

    sleep ${params.sleepTimeEnd}

    """
}

// run cellranger count

process cellranger_count {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus params.countJobCpus
    executor 'slurm'
    queue 'htc'
    time params.countJobLength
    memory params.countJobMemory

    tag "Running cellranger count on the samples"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    set sampleID, file(read) from sampleNames
    file organism from mkref_ch

    output:
    file("cellranger_count_${sampleID}") into count_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    ${cellrangerPath}/cellranger count --id=cellranger_count_${sampleID} --fastqs=${params.data} --sample=${sampleID} --transcriptome=${organism} --localcores 12 --localmem 80

    sleep ${params.sleepTimeEnd}
    """
}

