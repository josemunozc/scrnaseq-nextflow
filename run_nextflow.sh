#!/bin/bash

#SBATCH -p htc
#SBATCH --mem=80GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH -t 08:00:00 
#SBATCH -o slurm/%J.out
#SBATCH -e slurm/%J.err
#SBATCH --job-name=NEXT
#SBATCH --account=scw1448

myDir=$(pwd)

module load java 
module load nextflow/24.10.4
module load cellranger/7.2.0

#/scratch/c.mcbsd1/nf-core/scrnaseq/bin/nextflow-24.10.2/nextflow run nf-core/scrnaseq -profile singularity --input samplesheet.csv --cellranger_index /scratch/c.mcbsd1/nf-core/scrnaseq/resources/refdata-gex-mm10-2020-A --protocol 10XV2 --aligner cellranger --outdir ../output

nextflow run ${myDir}/main_dsl2_new2_final.nf --aggr true -with-trace


# to run
# sbatch run_nextflow.sh


