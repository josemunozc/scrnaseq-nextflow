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


Channel.fromPath("${params.dataDir}/targets.csv")
       .splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .into { sampleNames }


// concatenate and rename the input fastq files

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
    time params.concatJobLength
    memory params.concatJobMemory

    module params.ravenVersion
    module params.python3Version

    tag "Concatenating and renaming input fastq files"
    publishDir path:{params.outputDir},mode: 'copy'

    tag "Making gff from gtf for dexseq"
    publishDir path:{params.resourcesDir},mode: 'copy'

    input:
    file genome from genome_ch
    file gtf from gtf_ch

    output:
    file("${params.organism}") into mkref_ch

    script:
    """

    cellranger mkref --genome=${params.organism} --fasta=${genome} --genes=${gtf}

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

    cellranger count --id=cellranger_count_${sampleID} --fastqs=${params.data} --sample=${sampleID} --transcriptome=${organism} --localcores 24 --localmem 80

    sleep ${params.sleepTimeEnd}
    """
}



process trimgalore {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus params.trimgaloreJobCpus
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.trimgaloreJobLength
    memory params.trimgaloreJobMemory

    module params.ravenVersion
    module params.trimgaloreVersion

    tag "Running trim galore on the samples"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from concatfastq1_ch

    output:
    file("${sampleID}") into trimmedfastq1_ch
    file("${sampleID}") into trimmedfastq2_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    trim_galore --cores ${params.trimgaloreJobCpus} -o ${sampleID} --paired ${sampleID}/${sampleID}_1.fastq.gz ${sampleID}/${sampleID}_2.fastq.gz

    mv ${sampleID}/${sampleID}_1_val_1.fq.gz ${sampleID}/${sampleID}.trimmed_1.fastq.gz

    mv ${sampleID}/${sampleID}_2_val_2.fq.gz ${sampleID}/${sampleID}.trimmed_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """
}

// run fastqc on raw fastq

process raw_fastqc {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.fastqcJobLength
    memory params.fastqcJobMemory

    module params.fastqcVersion

    tag "Running fastQC on raw fastq"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from concatfastq2_ch

    output:
    file("${sampleID}") into fastqonraw_ch

    script:
    """ 
    sleep ${params.sleepTimeStart}

    fastqc -o ${sampleID} ${sampleID}/${sampleID}_1.fastq.gz

    fastqc -o ${sampleID} ${sampleID}/${sampleID}_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """

}

// run fastqc on trimmed fastq

process trimmed_fastqc {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.fastqcJobLength
    memory params.fastqcJobMemory

    module params.fastqcVersion

    tag "Running fastQC on trimmed fastq"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from trimmedfastq1_ch

    output:
    file("${sampleID}") into fastqontrimmed_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_1.fastq.gz

    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """

}

// run multiqc on all fastqc results

process multiqc {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    queue 'c_compute_wgp' 
    time params.multiqcJobLength
    memory params.multiqcJobMemory

    module params.ravenVersion
    module params.multiqcVersion

    tag "Running multiqc on fastqc reports"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy from fastqontrimmed_ch.collect()

    output:
    file("multiQC.html") into multiqc1_ch
    file("multiQC.html") into multiqc2_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir trimmed_fastqc

    multiqc \$linkDir -n multiQC

    rm -r \$linkDir

    sleep ${params.sleepTimeEnd}
    """
}

// get sequence length for STAR parameter

process sequencelength {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    executor 'slurm'
    queue 'c_compute_wgp' 
    time '10m'

    tag "Calculating sequence length for STAR genome index"
    publishDir path:{params.resourcesDir},mode: 'copy'

    output:
    file("seqlen.txt") into sequencelength_ch

    script:
    """
    bash ${params.projectDir}/${params.srcDir}/calculate_fastq_sequence_length.sh ${params.dataDir} > seqlen.txt
    """
}

// copy genome assembly file into resources

process copy_genome {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp' 
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

// copy GTF file into resources

process copy_gtf {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp' 
    time params.copyGTFJobLength

    tag "Copying gtf file into resources/"
    publishDir path:{params.resourcesDir},mode: 'copy'

    output:
    file("*") into gtf1_ch
    file("*") into gtf2_ch
    file("*") into gtf3_ch
    file("*") into gtf4_ch
    file("*") into gtf5_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    cp ${params.gtf} .

    sleep ${params.sleepTimeEnd}
    """
}

// create gff from gtf for dexseq

process make_gff {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.makeGFFJobLength
    memory params.makeGFFJobMemory

    module params.ravenVersion
    module params.htseqVersion
    module params.python2Version
    module params.dexseqVersion

    tag "Making gff from gtf for dexseq"
    publishDir path:{params.resourcesDir},mode: 'copy'

    input:
    file gtf from gtf1_ch

    output:
    file("genome.gff") into gff1_ch
    file("genome.gff") into gff2_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    python2.7 ${params.dexseq_prepare_annotationExecutable} ${gtf} genome.gff

    sleep ${params.sleepTimeEnd}
    """
}

// create STAR genome index

process genome_index {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus params.genomeIndexJobCpus
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.genomeIndexJobLength
    memory params.genomeIndexJobMemory

    module params.ravenVersion
    module params.starVersion

    tag "Indexing genome using STAR"
    publishDir path:{params.resourcesDir},mode: 'copy'

    input:
    file gtf from gtf2_ch
    file genome from genome_ch
    file seqlen from sequencelength_ch

    output:
    file "genome_index" into genomeindex_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    mkdir genome_index

    STAR --runThreadN ${params.genomeIndexJobCpus} --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf} --sjdbOverhang `cat ${seqlen}`

    sleep ${params.sleepTimeEnd}
    """
}

// run STAR mapping

process star_mapping {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus params.starMappingJobCpus
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.starMappingJobLength
    memory params.starMappingJobMemory

    module params.ravenVersion
    module params.starVersion

    tag "Mapping trimmed fastq with STAR"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from trimmedfastq2_ch
    file genomeIndex from genomeindex_ch

    output:
    file("${sampleID}") into bam_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    STAR --readFilesCommand zcat --runThreadN ${params.starMappingJobCpus} --runMode alignReads --outSAMunmapped Within KeepPairs --outMultimapperOrder Random --outSAMmultNmax 1 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sampleID}/${sampleID}.randomonemap. --genomeDir ${genomeIndex} --readFilesIn ${sampleID}/${sampleID}.trimmed_1.fastq.gz ${sampleID}/${sampleID}.trimmed_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """
}

// mark duplicate reads in STAR bam

process markduplicates {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 4
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.markDuplicatesJobLength
    memory params.markDuplicatesJobMemory

    module params.ravenVersion
    module params.picardVersion
    module params.samtoolsVersion

    tag "Running markduplicates sorted bams using picard"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from bam_ch

    output:
    file("${sampleID}") into markdupbam1_ch
    file("${sampleID}") into markdupbam2_ch
    file("${sampleID}") into markdupbam3_ch
    file("${sampleID}") into markdupbam4_ch
    file("${sampleID}") into markdupbam5_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    java -jar -Xmx40G ${params.picardExecutable} MarkDuplicates I=${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam O=${sampleID}/${sampleID}.markdup.bam M=${sampleID}/${sampleID}.markdup.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

    java -jar -Xmx40G ${params.picardExecutable} MarkDuplicates I=${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam O=${sampleID}/${sampleID}.rmdup.bam M=${sampleID}/${sampleID}.rmdup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

    samtools index ${sampleID}/${sampleID}.markdup.bam

    samtools index ${sampleID}/${sampleID}.rmdup.bam

    samtools sort -n ${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam -o ${sampleID}/${sampleID}.markdup.namesorted.bam

    samtools sort -n ${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam -o ${sampleID}/${sampleID}.rmdup.namesorted.bam

    sleep ${params.sleepTimeEnd}
    """
}

// run bamtools stats on marked bams

process bamtools {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    executor 'slurm'
    queue 'c_compute_wgp'
    time params.bamtoolsJobLength
    memory params.bamtoolsJobMemory

    module params.ravenVersion
    module params.bamtoolsVersion

    tag "Running bamtools on markdup bams"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from markdupbam1_ch

    output:
    file("${sampleID}") into bamtools_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    bamtools stats -in ${sampleID}/${sampleID}.markdup.bam > ${sampleID}/${sampleID}.markdup.stats.txt

    sleep ${params.sleepTimeEnd}
    """

}

// run picard stats

process collect_insert_size_metrics {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'

    time params.collectInsertSizeMetricsJobLength
    memory params.collectInsertSizeMetricsJobMemory

    module params.picardVersion
    module params.ravenVersion
    module params.picardRVersion

    tag "Running picard collect insert size metrics on markdup bams"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from bamtools_ch

    output:
    file("${sampleID}") into picardmetrics_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    java -jar -Xmx40G ${params.picardExecutable} CollectInsertSizeMetrics I=${sampleID}/${sampleID}.markdup.bam H=${sampleID}/${sampleID}.markdup.picardstats.out O=${sampleID}/${sampleID}.markdup.picardstats.txt

    sleep ${params.sleepTimeEnd}
    """
}

// make bamtools report

process run_bamtools_report {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.bamtoolsReportJobLength

    tag "Collating results from bamtools"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy1 from picardmetrics_ch.collect()

    output:
    file("*.txt") into bamtoolsreport_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir raw_fastqc
    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir bamtools
    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir collect_insert_size_metrics 

    ${params.projectDir}/${params.srcDir}/collate_bamtools_and_picard_output.pl ${params.dataDir}/targets.csv ${params.projectDir} ${params.outputDir} \$linkDir

    rm -r \$linkDir

    sleep ${params.sleepTimeEnd}
    """
}

// count reads for all genes and transcripts using featureCounts

process featurecounts {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus params.featureCountsJobCpus 
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.featureCountsJobLength
    memory params.featureCountsJobMemory

    module params.featureCountsVersion

    tag "Running featurecounts for all genes and transcripts"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from markdupbam2_ch
    file gtf from gtf3_ch

    output:
    file("${sampleID}") into featurecounts_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g gene_id -a ${gtf} -o ${sampleID}/${sampleID}.markdup.genecount.txt ${sampleID}/${sampleID}.markdup.bam

    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g gene_id -a ${gtf} -o ${sampleID}/${sampleID}.rmdup.genecount.txt ${sampleID}/${sampleID}.rmdup.bam

    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g transcript_id -a ${gtf} -o ${sampleID}/${sampleID}.markdup.transcriptcount.txt ${sampleID}/${sampleID}.markdup.bam

    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g transcript_id -a ${gtf} -o ${sampleID}/${sampleID}.rmdup.transcriptcount.txt ${sampleID}/${sampleID}.rmdup.bam

    sleep ${params.sleepTimeEnd}
    """
}

// create counts file for dexseq input

process dexseq_count {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }

    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.dexseqCountJobLength
    memory params.dexseqCountJobMemory

    module params.ravenVersion
    module params.htseqVersion
    module params.python2Version
    module params.dexseqVersion
    module params.samtoolsVersion

    tag "Running creating counts for dexseq input"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file sampleID from markdupbam3_ch
    file gff from gff1_ch

    output:
    file("${sampleID}") into dexseqcount1_ch
    file("${sampleID}") into dexseqcount2_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    samtools view -F 256 ${sampleID}/${sampleID}.markdup.namesorted.bam | python2.7 ${params.dexseq_countExecutable} -p yes ${gff} - ${sampleID}/${sampleID}.markdup.dexseqcount.txt

    sleep ${params.sleepTimeEnd}
    """
}

// make count report

process run_featurecounts_report {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.featureCountsReportJobLength
    memory params.featureCountsReportJobMemory

    tag "Collating results from featureCounts"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy from featurecounts_ch.collect()
    file gtf from gtf4_ch

    output:
    file("all*.txt") into featurecountsreport_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir featurecounts

    ${params.projectDir}/${params.srcDir}/collate_featurecounts_output.pl ${params.dataDir}/targets.csv ${params.projectDir} ${params.outputDir} ${params.resourcesDir} ${gtf} \$linkDir

    rm -r \$linkDir

    sleep ${params.sleepTimeEnd}
    """
}

// run deseq2 normalising over all data

process deseq2_allnorm {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 4
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.deseq2JobLength

    module params.raven
    module params.rVersion

    tag "Running deseq2 all norm"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy from featurecountsreport_ch.collect()

    output:
    file("deseq2-allnorm") into deseq2allnorm_ch

    script:
    """
    mkdir -p deseq2-allnorm/output
    mkdir -p deseq2-allnorm/qc

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir run_featurecounts_report

    Rscript ${params.projectDir}/${params.srcDir}/deseq2-allnorm.R \$linkDir ${params.dataDir}/targets.csv

    rm -r \$linkDir

    """
}

// run deseq2 normalising within the comparison pair

process deseq2_pairnorm {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 4
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.deseq2JobLength

    module params.raven
    module params.rVersion

    tag "Running deseq2 pair norm"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy from deseq2allnorm_ch 

    output:
    file("deseq2-pairnorm") into deseq2pairnorm_ch

    script:
    """
    mkdir -p deseq2-pairnorm/output
    mkdir -p deseq2-pairnorm/qc

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir run_featurecounts_report

    Rscript ${params.projectDir}/${params.srcDir}/deseq2-pairnorm.R \$linkDir ${params.dataDir}/targets.csv

    rm -r \$linkDir

    """
}

// run dexseq 

process dexseq {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    clusterOptions params.dexseqClusterOptions
    executor 'slurm'

    module params.ravenVersion
    module params.rVersion
    module params.dexseqVersion

    tag "Running dexseq"
    publishDir path:{params.outputDir},mode: 'copy'

    input:
    file dummy1 from dexseqcount1_ch.collect()
    file dummy2 from gff2_ch
    file dummy3 from deseq2pairnorm_ch

    output:
    file("dexseq") into dexseq_ch

    script:
    """
    mkdir -p dexseq/output

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir dexseq_count

    Rscript ${params.projectDir}/${params.srcDir}/dexseq.R genome.gff ${params.projectDir}/${params.resourcesDir} ${params.dataDir}/targets.csv \$linkDir 6 2> /dev/null

    rm -r \$linkDir

    """
}

// make release directory

process make_release {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1

    cpus 1
    executor 'slurm'
    queue 'c_compute_wgp'
    time params.makeReleaseJobLength

    tag "Making release directory"
    publishDir path:{params.releaseDir},mode: 'copy'

    input:
    file dummy01 from bamtoolsreport_ch.collect()
    file dummy02 from featurecountsreport_ch.collect()
    file dummy03 from markdupbam4_ch.collect()
    file dummy04 from multiqc1_ch
    file dummy05 from dexseq_ch

    output:
    file("*") into release_ch

    script:
    """

    ${params.projectDir}/${params.srcDir}/make_release.pl ${params.projectDir} ${params.analysisID} targets.csv ${params.dataDir} ${params.outputDir} ${params.resourcesDir} ${params.gtfName}


    """
}





