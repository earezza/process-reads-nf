#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-jdilwort
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=nextflow-example
#SBATCH --output=%x-%j.out

module load nextflow
module load apptainer

export NXF_DISABLE_CHECK_LATEST=true

cd ~/scratch/nextflows/

nextflow run process_reads.nf -c process_reads.config -profile cluster \
	-with-dag -with-report -with-timeline -w ~/scratch/nextflows/work-cnt/ \
	--target CUTnTag/ \
	--reads_type paired \
	--assembly mm10 \
	--length 100 \
	--assay cutntag \
	--multiqc_config multiqc_config.yaml \
	--genome_index ~/projects/def-jdilwort/Reference_Files/bowtie2/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
