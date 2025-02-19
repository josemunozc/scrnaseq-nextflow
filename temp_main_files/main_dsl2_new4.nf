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
params.index = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/index"

workflow {

    read_pairs_ch = Channel.fromPath("${params.dataDir}/samplesheet_1kmouse.csv").splitCsv(header:true)
       .map{ row-> tuple(row.suppliedID) }
       .view()
    CRCOUNT(read_pairs_ch)
}
// copy genome files

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
    tuple val(read_pairs_ch)

    output:
    path "cellranger_count_${sampleID}"

    script:
    """
    sleep ${params.sleepTimeStart}

    ${params.cellrangerPath}/cellranger count --id=cellranger_count_${read_pairs_ch} --fastqs=${params.data} --sample=${read_pairs_ch} --transcriptome=${params.index} --localcores 12 --localmem 80

    sleep ${params.sleepTimeEnd}
    """
}

