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

params.outdir = "results"
params.resourcesDir = "resources"
params.genome = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.cdna.all.fa"
params.gtf = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.113.gtf"
params.organism = "Mus_musculus"
params.dataDir = "/scratch/c.mcbsd1/scrnaseq-nextflow"
params.cellrangerPath = "/scratch/c.mcbsd1/cellranger/cellranger-7.1.0"

workflow {

    read_pairs_ch = Channel.fromPath("${params.dataDir}/samplesheet_1kmouse.csv").splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .view()
    MKREF(params.cellrangerPath, params.genome, params.gtf)
    CRCOUNT(read_pairs_ch, MKREF.out)
}
// copy genome files

process MKREF {

    maxRetries params.retries
    maxErrors -1

    cpus 1  
    executor 'slurm'
    queue 'htc'
    time params.mkrefJobLength
    memory params.mkrefMemory

    tag "Making reference files"
    publishDir params.resourcesDir

    input:
    path genome
    path gtf
    path cellrangerPath

    output:
    path 'index'

    script:
    """

    ${cellrangerPath}/cellranger mkref --genome=index --fasta=${genome} --genes=${gtf}

    sleep ${params.sleepTimeEnd}

    """
}

// run cellranger count

process CRCOUNT {

    maxRetries params.retries
    maxErrors -1

    cpus params.countJobCpus
    executor 'slurm'
    queue 'htc'
    time params.countJobLength
    memory params.countJobMemory

    tag "Running cellranger count on the samples"
    publishDir params.outdir

    input:
    tuple val(sampleID)
    path index

    output:
    path "cellranger_count_${sampleID}"

    script:
    """
    sleep ${params.sleepTimeStart}

    ${cellrangerPath}/cellranger count --id=cellranger_count_${sampleID} --fastqs=${params.data} --sample=${sampleID} --transcriptome=${organism} --localcores 12 --localmem 80

    sleep ${params.sleepTimeEnd}
    """
}

