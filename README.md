Nextflow implementation for processing Illumina short reads from high-throughput assays such as CUT&Tag and RNA-Seq.  
Apptainer containing all software required to run each process can be built using the dilworthlab.def file:  
    sudo apptainer build --mksquashfs-args="-comp gzip" --force dilworthlab.sif dilworthlab.def    

Modify process_reads.config to point full path to your dilworthlab.sif file (container = ...).  

Then can run, e.g.  
    nextflow run process_reads.nf -c process_reads.config \
        -with-dag -with-report -with-timeline \
        --target Sample/ \
        --reads_type paired \
        --assembly mm10 \
        --length 150 \
        --normalize_by BPM \
        --assay rnaseq \
        --multiqc_config multiqc_config.yaml \
        --genome_index hisat2/mm10/genome 
