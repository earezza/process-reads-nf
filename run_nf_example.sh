#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-jdilwort
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --array=0-3
#SBATCH --job-name=NF-EXAMPLE
#SBATCH --output=%x-%j.out

module load nextflow
module load apptainer

export NXF_DISABLE_CHECK_LATEST=true

cd ~/scratch/

samples=(Control-1 Control-2 KO-1 KO-2)

nextflow run process_reads.nf -c process_reads.config \
	-with-dag -with-report -with-timeline \
	--target ${samples[SLURM_ARRAY_TASK_ID]}/ \
	--reads_type paired \
	--assembly mm10 \
	--length 100 \
	--assay cutntag \
	--multiqc_config multiqc_config.yaml \
	--genome_index Reference_Files/bowtie2/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
        --ecoli_index Reference_Files/Spikein_indices/EcoliK12_index/EcoliK12Index/EcoliK12
