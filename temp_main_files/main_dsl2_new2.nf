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

params.workingDir = "/scratch/c.mcbsd1/scrnaseq-nextflow"
params.outdir = "results"
params.resourcesDir = "resources"
params.genome = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
params.gtf = "/scratch/c.mcbsd1/scrnaseq-nextflow/resources/Mus_musculus.GRCm39.113.gtf"
params.organism = "Mus_musculus"
params.data = "/scratch/c.mcbsd1/scrnaseq-nextflow/data"
params.cellrangerPath = "/scratch/c.mcbsd1/cellranger/cellranger-7.1.0"
params.aggr = false
params.res = "/scratch/c.mcbsd1/scrnaseq-nextflow/results"

workflow {

    read_pairs_ch = Channel.fromPath("${params.workingDir}/samplesheet_1kmouse.csv").splitCsv(header:true)
       .map{ row-> tuple(row.suppliedID) }
       .view()
    MKREF()
    CRCOUNT(read_pairs_ch, MKREF.out)
    if (params.aggr) {
      GENAGGR(read_pairs_ch,CRCOUNT.out)
      AGGRCOUNT(read_pairs_ch,GENAGGR.out)
    }
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

    output:
    path 'index', emit: index

    script:
    """

    ${params.cellrangerPath}/cellranger mkref --genome=index --fasta=${params.genome} --genes=${params.gtf}

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
    tuple val(read_pairs_ch)
    path index

    output:
    path "cellranger_count_${read_pairs_ch}"

    script:
    """
    sleep ${params.sleepTimeStart}

    ${params.cellrangerPath}/cellranger count --id=cellranger_count_${read_pairs_ch} --fastqs=${params.data} --sample=${read_pairs_ch} --transcriptome=${index} --localcores 12 --localmem 80

    sleep ${params.sleepTimeEnd}
    """
}

process GENAGGR {
    tag "Generate aggr.csv file"
    publishDir params.outdir

    input:
    tuple val(read_pairs_ch)
    path "cellranger_count_${read_pairs_ch}"

    output:
    path "aggr.csv"

    script:
    """
    echo "sample_id,molecule_h5" > aggr.csv
    echo "${read_pairs_ch},${params.res}/cellranger_count_${read_pairs_ch}/outs/molecule_info.h5" >> aggr.csv
    """
}

process AGGRCOUNT {
    tag "Aggregate cellranger count results"
    publishDir params.outdir

    input:
    tuple val(read_pairs_ch)    
    path "aggr.csv"

    output:
    path "aggregation"

    script:
    """
    ${params.cellrangerPath}/cellranger aggr --id=aggregation --csv=aggr.csv
    """
}
