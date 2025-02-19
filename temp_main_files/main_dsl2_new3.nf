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
params.genome = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
params.gtf = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.113.gtf"
params.organism = "Mus_musculus"
params.dataDir = "/scratch/c.mcbsd1/scrnaseq-nextflow"
params.cellrangerPath = "/scratch/c.mcbsd1/cellranger/cellranger-7.1.0"

workflow {

    read_pairs_ch = Channel.fromPath("${params.dataDir}/samplesheet_1kmouse.csv").splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .view()
    MKREF()
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

    output:
    path 'index', emit: index

    script:
    """

    ${params.cellrangerPath}/cellranger mkref --genome=index --fasta=${params.genome} --genes=${params.gtf}

    sleep ${params.sleepTimeEnd}

    """
}

