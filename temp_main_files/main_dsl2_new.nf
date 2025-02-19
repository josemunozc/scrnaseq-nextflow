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
params.genome = "Mus_musculus.GRCm39.cdna.all.fa"
params.gtf = "Mus_musculus.GRCm39.113.gtf"
params.organism = "Mus_musculus"
params.dataDir = "/scratch/c.mcbsd1/scrnaseq-nextflow"

workflow {

    read_pairs_ch = Channel.fromPath("${params.dataDir}/samplesheet_1kmouse.csv")
       .splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .into { sampleNames }
    MKREF(params.genome, params.gtf, COPYGENOME.out, COPYGTF.out)
    CRCOUNT(read_pairs_ch, MKREF.out)
}
// copy genome files

process MKREF {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
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
    path organism

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

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
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
    tuple sampleID, file(read) from sampleNames
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

