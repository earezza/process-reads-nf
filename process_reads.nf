#!/usr/bin/env nextflow

def display_start() {
    log.info """\
        ===================================================
            P R O C E S S   R E A D S   P I P E L I N E
        ===================================================
        DESCRIPTION:
        Process Illumina raw reads (.fastq.gz) from CUT&Tag, 
        ChIPseq, ATACseq, RNAseq, or other high-throughput 
        next-generation genomic/transcriptomic sequences.
        OUTPUTS:
        - QC reports (.html)
        - Alignment (.bam)
        - Coverage (.bw)
        - Peaks (.bed, .narrowPeak, .broadPeak)
        ---------------------------------------------------
        AUTHOR: 
        Eric Arezza
        earezza17@gmail.com
        ---------------------------------------------------
        VERSION:
        0.0.1
        ===================================================
        target          : ${params.target}
        control         : ${params.control}
        reads_type      : ${params.reads_type}
        length          : ${params.length}
        assay           : ${params.assay}
        assembly        : ${params.assembly}
    """.stripIndent(true)
}

/*
    ---------------------------------------------------------
        PREP COMMANDS AND VARIABLES BASED ON INPUT PARAMS
    ---------------------------------------------------------
*/
def build_params(){

    // Reference index file for spike-in
    /*
    if (params.spikein_type == 'Amp'){
        params.spike_index = params.spikein_index_amp
    } else if (params.spikein_type == 'Bacteria'){
        params.spike_index = params.spikein_index_ecoli
    } else {
        throw new Exception ("Spike-in index not available")
    }
    */
}

/*
    -------------------------------------------------
        VERIFY AND VALIDATE READS FILES FOR INPUT
    -------------------------------------------------
*/
// To load reads files
def get_reads_ch(reads_dir){
    if ("${reads_dir}" != ''){
        assert file("${reads_dir}").exists() : "Cannot find target reads directory ${reads_dir}"
        if (params.reads_type == 'paired'){
            assert files("${reads_dir}/*_R{1,2}.fastq.gz").size() > 0 : "Expected paired reads...check files in ${reads_dir}"
            assert files("${reads_dir}/*_R{1,2}.fastq.gz").size() %2 == 0 : "Expected paired reads...check files in ${reads_dir}"
            def reads = Channel.fromFilePairs("${reads_dir}/*_R{1,2}.fastq.gz", type: 'file')
                .map { items -> tuple(file("${reads_dir}"), items[1]) }
                //.groupTuple()
            return reads
        } else {
            assert files("${reads_dir}/*.fastq.gz").size() > 0 : "Expected reads...check files in ${reads_dir}"
            def reads = Channel.fromPath("${reads_dir}/*.fastq.gz", type: 'file')
                .map { items -> tuple(file("${reads_dir}"), items) }
                .groupTuple()
            return reads
        }
    } else {
        def reads = Channel.empty()
        return reads
    }
}

// To load supplemental files
/*
def get_file_ch( file_path ) {
    if ("${file_path}" != ''){
        if (file( "${file_path}" ).exists()){
            def file_ch = Channel.fromPath( "${file_path}" , type: 'file', checkIfExists: true)
            return file_ch
        } else {
            def file_ch = Channel.empty()
            return file_ch
        }
    } else {
        def file_ch = Channel.empty()
        return file_ch
    }
}
*/
include { MD5SUMCHECK; MULTIQC; FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED } from './modules.nf'
include { CUTADAPT_PAIRED; CUTADAPT_SINGLE } from './modules.nf'
include { HISAT2_PAIRED; HISAT2_SINGLE; BOWTIE2_PAIRED; BOWTIE2_SINGLE } from './modules.nf'
include { MERGE_BAMS; SORT_BAM; QFILTER_BAM; DEDUPLICATE_BAM; INDEX_BAM} from './modules.nf'
include { BIGWIG_COVERAGE; BIGWIG_BAMCOMPARE; BIGWIG_COVERAGE_STRANDED as BIGWIG_COVERAGE_STRANDED_FORWARD; BIGWIG_COVERAGE_STRANDED as BIGWIG_COVERAGE_STRANDED_REVERSE; BEDGRAPH_COVERAGE } from './modules.nf'
include { MACS_NARROWPEAKS_CONTROL; MACS_BROADPEAKS_CONTROL; MACS_NARROWPEAKS_NO_CONTROL; MACS_BROADPEAKS_NO_CONTROL; } from './modules.nf'
include { SEACR_PEAKS_CONTROL; SEACR_PEAKS_NO_CONTROL } from './modules.nf'
include { GOPEAKS_PEAKS_CONTROL; GOPEAKS_PEAKS_NO_CONTROL } from './modules.nf'

// Define workflow
workflow {

    display_start()

    // Read lengths genome mapping based on https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    def mm10_egs = [50: '2701495761', 75: '2747877777', 100: '2805636331', 150: '2862010578', 200: '2887553303']
    def hg38_egs = [50: '2308125349', 75: '2407883318', 100: '2467481108', 150: '2494787188', 200: '2520869189']
    def rn6_egs = [50: '2375372135', 75: '2440746491', 100: '2480029900', 150: '2477334634', 200: '2478552171']
    if (params.assembly == 'mm10'){
        effective_genome_size = Channel.value( mm10_egs[ params.length ] )
    } else if (params.assembly == 'hg38'){
        effective_genome_size = Channel.value( hg38_egs[ params.length ] )
    } else if (params.assembly == 'rn6'){
        effective_genome_size = Channel.value( rn6_egs[ params.length ] )
    } else {
        throw new Exception ("Assembly not available to estimate effective_genome_size")
    }

    // Load .fastq.gz files
    target_reads = get_reads_ch( params.target )
    control_reads = get_reads_ch( params.control )

    // Get .md5 files
    target_reads
        .flatMap{ id, files ->
            files.collect { item -> tuple(id, item, file("${item}.md5")) }
            }
        .set{ target_md5 }
    control_reads
        .flatMap{ id, files ->
            files.collect { item -> tuple(id, item, file("${item}.md5")) }
            }
        .set{ control_md5 }

    // Collect replicates of same sample into separate processing channels
    if (params.reads_type == 'single'){
        target_reads
            .flatMap{ sample, reads -> reads.collect{ read -> tuple( sample, read)  } }
            .set{ target_reads }
        control_reads
            .flatMap{ sample, reads -> reads.collect{ read -> tuple( sample, read ) } }
            .set{ control_reads }
    }/* else {
        Channel.fromFilePairs("${params.target}/*_R{1,2}.fastq.gz", type: 'file')
            .map { items -> tuple(file("${params.target}"), items[1]) }
            .groupTuple()
            .set{ target_reads }
        Channel.fromFilePairs("${params.control}/*_R{1,2}.fastq.gz", type: 'file')
            .map { items -> tuple(file("${params.control}"), items[1]) }
            .groupTuple()
            .set{ control_reads }
    }*/
    
    // Concat reads list
    target_reads
        .concat( control_reads )
        .set { raw_reads }

    
    // Check file integrity
    MD5SUMCHECK( target_md5.concat( control_md5 ) )

    // QC on raw reads
    FASTQC_RAW( raw_reads )

    // Process raw reads
    if (params.reads_type == 'paired'){
        // Trim
        CUTADAPT_PAIRED( raw_reads )
        // QC on trimmed reads
        FASTQC_TRIMMED( CUTADAPT_PAIRED.out.trimmed_reads )
        // Map reads to reference assembly and sort alignment
        if (params.assay == 'rnaseq') {
            HISAT2_PAIRED( CUTADAPT_PAIRED.out.trimmed_reads )
            SORT_BAM( HISAT2_PAIRED.out.mapped_reads )
        } else {
            BOWTIE2_PAIRED( CUTADAPT_PAIRED.out.trimmed_reads )
            SORT_BAM( BOWTIE2_PAIRED.out.mapped_reads )
        }
    } else {
        // Trim
        CUTADAPT_SINGLE( raw_reads )
        // QC on trimmed reads
        FASTQC_TRIMMED( CUTADAPT_SINGLE.out.trimmed_reads )
        // Map reads to reference assembly and sort alignment
        if (params.assay == 'rnaseq') {
            HISAT2_SINGLE( CUTADAPT_SINGLE.out.trimmed_reads )
            SORT_BAM( HISAT2_SINGLE.out.mapped_reads )
        } else {
            BOWTIE2_SINGLE( CUTADAPT_SINGLE.out.trimmed_reads )
            SORT_BAM( BOWTIE2_SINGLE.out.mapped_reads )
        }
    }

    // Filter out poor alignments
    QFILTER_BAM( SORT_BAM.out.bam_sorted )

    // Remove duplicates
    DEDUPLICATE_BAM( QFILTER_BAM.out.bam_qfiltered )

    // Index final alignments for filtered and deduplicated bams
    INDEX_BAM( QFILTER_BAM.out.bam_qfiltered.concat(DEDUPLICATE_BAM.out.bam_deduplicated) )

    // Separate target and control
    if (params.assay == 'rnaseq') {
        // Use alignments containing duplicates for RNAseq
        INDEX_BAM.out.bam_indexed
        .branch {
            target: it[0].baseName == params.target.replace('/', '') && !it[1].name.contains('NoDups')
            control: it[0].baseName == params.control.replace('/', '') && !it[1].name.contains('NoDups')
        }.set { bams }
    } else {
        // Use alignments with duplicates removed for other assays
        INDEX_BAM.out.bam_indexed
        .branch {
            target: it[0].baseName == params.target.replace('/', '') && it[1].name.contains('NoDups')
            control: it[0].baseName == params.control.replace('/', '') && it[1].name.contains('NoDups')
        }.set { bams }
    }

    // TODO: Merge replicate alignments into a representative for the sample (typically used when multiple single-end reads in SAMPLE/ directory)
    if (params.merge) {
        bams.target.groupTuple().concat(bams.control.groupTuple()).view()
        //error(" Exiting script ")
        //MERGE_BAMS( bams.target.groupTuple().concat(bams.control.groupTuple()) )

        /*
        MERGE_BAMS.out.bam_merged
        .branch {
            target_merged: it[0].baseName == params.target.replace('/', '')
            control_merged: it[0].baseName == params.control.replace('/', '')
        }.set { bams_merged }
        */
    }

    // Run MultiQC report
    multiqc_config  = (params.multiqc_config != '') ? Channel.value( file(params.multiqc_config) ) : Channel.value( file('NO_MULTIQC_CONFIG_FILE') )
    MULTIQC( bams.target.concat(bams.control), multiqc_config )
    
    // Get genome coverage bigwigs
    blacklist_regions  = (params.blacklist != '') ? Channel.value( file(params.blacklist) ) : Channel.value( file('NO_BLACKLIST_FILE') )
    BIGWIG_COVERAGE( bams.target.concat( bams.control ), blacklist_regions, effective_genome_size )
    if (params.control != ''){
        BIGWIG_BAMCOMPARE( bams.target, bams.control, blacklist_regions, effective_genome_size )
    }
    
    // Peak calling or stranded bigwigs if RNAseq
    if (params.assay != 'rnaseq'){
        // bedgraphs for seacr
        BEDGRAPH_COVERAGE( bams.target.concat(bams.control), blacklist_regions, effective_genome_size )
        BEDGRAPH_COVERAGE.out.bedgraph
            .branch {
                target: it[0].baseName == params.target.replace('/', '')
                control: it[0].baseName == params.control.replace('/', '')
        }.set { bedgraphs }
        // Call peaks
        // Peaks independently
        MACS_NARROWPEAKS_NO_CONTROL( bams.target.concat(bams.control), effective_genome_size )
        MACS_BROADPEAKS_NO_CONTROL( bams.target.concat(bams.control), effective_genome_size )
        SEACR_PEAKS_NO_CONTROL( bedgraphs.target.concat(bedgraphs.control) )
        GOPEAKS_PEAKS_NO_CONTROL( bams.target.concat(bams.control) )
        if (params.control != ''){
            // Peaks with control
            MACS_NARROWPEAKS_CONTROL( bams.target, bams.control, effective_genome_size )
            MACS_BROADPEAKS_CONTROL( bams.target, bams.control, effective_genome_size )
            SEACR_PEAKS_CONTROL( bedgraphs.target, bedgraphs.control )
            GOPEAKS_PEAKS_CONTROL( bams.target, bams.control )
        }
    } else {
        if (params.stranded){
            BIGWIG_COVERAGE_STRANDED_FORWARD ( INDEX_BAM.out.bam_indexed, blacklist_regions, effective_genome_size, 'forward' )
            BIGWIG_COVERAGE_STRANDED_REVERSE ( INDEX_BAM.out.bam_indexed, blacklist_regions, effective_genome_size, 'reverse' )
        }
    }
    
    // TODO: handle spike-ins if used
    // Map to spike-in index
    // Calculate scaling/normalization factor
    // Get spike-in-normalized coverage
    
}

