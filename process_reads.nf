#!/usr/bin/env nextflow
/*
 * ===================================================
 *  Process Reads Pipeline - main workflow
 * ===================================================
 * Process Illumina raw reads (.fastq.gz) from CUT&Tag, ChIPseq, ATACseq,
 * RNAseq, or other high-throughput next-generation genomic/transcriptomic
 * sequencing.
 *
 * Run locally:  nextflow run test.nf -profile standard --target reads/ --genome_index /path/to/index
 * Run on SLURM: nextflow run test.nf -profile slurm    --target reads/ --genome_index /path/to/index
 */

include { MD5SUMCHECK; FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED; MULTIQC } from './modules.nf'
include { BOWTIE2 as ALIGN_MAIN; BOWTIE2 as ALIGN_ECOLI; BOWTIE2 as ALIGN_SPIKEIN } from './modules.nf'
include { HISAT2 as ALIGN_RNA_MAIN } from './modules.nf'
include { CUTADAPT; PREFILTER_QC; QFILTER_BAM; DEDUPLICATE_BAM } from './modules.nf'
include { GET_SCALE_FACTOR as ECOLI_SCALE_FACTOR_TARGET; GET_SCALE_FACTOR as ECOLI_SCALE_FACTOR_CONTROL; GET_SCALE_FACTOR as SPIKEIN_SCALE_FACTOR_TARGET; GET_SCALE_FACTOR as SPIKEIN_SCALE_FACTOR_CONTROL } from './modules.nf'
include { BIGWIG_COVERAGE as BIGWIGS_ECOLI_SCALE_FACTOR_TARGET; BIGWIG_COVERAGE as BIGWIGS_SPIKEIN_SCALE_FACTOR_TARGET; BIGWIG_COVERAGE as BIGWIGS_NO_SCALE_FACTOR_TARGET } from './modules.nf'
include { BIGWIG_COVERAGE as BIGWIGS_ECOLI_SCALE_FACTOR_CONTROL; BIGWIG_COVERAGE as BIGWIGS_SPIKEIN_SCALE_FACTOR_CONTROL; BIGWIG_COVERAGE as BIGWIGS_NO_SCALE_FACTOR_CONTROL } from './modules.nf'
include { BIGWIG_COVERAGE as BIGWIGS_TARGET_FORWARD; BIGWIG_COVERAGE as BIGWIGS_TARGET_REVERSE; BIGWIG_COVERAGE as BIGWIGS_CONTROL_FORWARD; BIGWIG_COVERAGE as BIGWIGS_CONTROL_REVERSE } from './modules.nf'
include { BIGWIG_COVERAGE as BIGWIGS_TARGET_UNSTRANDED; BIGWIG_COVERAGE as BIGWIGS_CONTROL_UNSTRANDED } from './modules.nf'
include { MERGE_BAMS; BIGWIG_COVERAGE as BIGWIGS_REPLICATES; BIGWIG_COVERAGE as BIGWIGS_REPLICATES_FORWARD; BIGWIG_COVERAGE as BIGWIGS_REPLICATES_REVERSE; BIGWIG_COVERAGE as BIGWIGS_REPLICATES_UNSTRANDED } from './modules.nf'
include { BIGWIG_BAMCOMPARE; BEDGRAPH_COVERAGE } from './modules.nf'
include { MACS_PEAKS; MACS_PEAKS_WITH_CONTROL; SEACR_PEAKS; SEACR_PEAKS_WITH_CONTROL; GOPEAKS_PEAKS; GOPEAKS_PEAKS_WITH_CONTROL } from './modules.nf'
include { NAME_SORT_BAM; GENRICH_PEAKS; GENRICH_PEAKS_WITH_CONTROL; FRIP } from './modules.nf'

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
        VERSION: ${workflow.manifest.version}
        NEXTFLOW: ${nextflow.version}
        PROFILE: ${workflow.profile}
        ===================================================
        target          : ${params.target}
        control         : ${params.control}
        reads_type      : ${params.reads_type}
        length          : ${params.length}
        assay           : ${params.assay}
        assembly        : ${params.assembly}
        genome_index    : ${params.genome_index}
        outdir          : ${params.outdir ?: params.target}
    """.stripIndent(true)
}

// Build a (sample_dir, reads) channel for a directory of fastq.gz files.
// One directory = one sample, so a directory must hold exactly one R1/R2 pair
// (paired-end) or one file (single-end) - every output is named after the
// directory. Exception: with --merge, a single-end directory may hold several
// files; each is processed as its own library (outputs named after the file)
// and the final BAMs are merged into the sample's BAM.
def get_reads_ch(reads_dir) {
    if (!reads_dir) {
        return channel.empty()
    }
    assert file(reads_dir).exists() : "Cannot find reads directory: ${reads_dir}"

    if (params.reads_type == 'paired') {
        def pairFiles = files("${reads_dir}/*_R{1,2}*.fastq.gz")
        assert pairFiles.size() > 0 : "Expected paired reads (R1/R2)...check files in ${reads_dir}"
        assert pairFiles.size() == 2 : "Expected exactly one R1/R2 pair in ${reads_dir} but found ${pairFiles.size()} files. Put each sample in its own directory (or concatenate lanes first)"
        return channel.fromFilePairs("${reads_dir}/*_R{1,2}*.fastq.gz", type: 'file')
            .map { _pairId, readPair -> tuple(file(reads_dir), readPair) }
    } else {
        def singleFiles = files("${reads_dir}/*.fastq.gz")
        assert singleFiles.size() > 0 : "Expected reads...check files in ${reads_dir}"
        assert params.merge || singleFiles.size() == 1 : "Found ${singleFiles.size()} single-end files in ${reads_dir}. Use --merge to process each as its own library and merge them into one sample, or put each sample in its own directory"
        if (params.merge) {
            // per-file outputs are named after the file, the merged ones after the directory
            def ids = singleFiles.collect { f -> f.name.replaceAll('(_trimmed)?\\.(fastq|fq)(\\.gz)?$', '') }
            assert ids.toSet().size() == ids.size() : "--merge: file names in ${reads_dir} must differ before .fastq.gz (found ${ids})"
            assert !ids.any { String id -> id.contains('.') } : "--merge: file names in ${reads_dir} may not contain '.' before .fastq.gz (found ${ids})"
            assert !(file(reads_dir).simpleName in ids) : "--merge: a file in ${reads_dir} is named like the directory itself ('${file(reads_dir).simpleName}'), so its outputs would clash with the merged ones. Rename the file"
        }
        return channel.fromPath("${reads_dir}/*.fastq.gz", type: 'file')
            .map { f -> tuple(file(reads_dir), f) }
    }
}

// Wrap an optional file path param, falling back to a sentinel file the
// downstream process scripts check for by name (e.g. NO_BLACKLIST_FILE).
// Returns a VALUE channel: these files are shared by every sample, and a
// one-item queue channel would let each process run only once (the first
// BAM to arrive would consume the single item and all later BAMs would wait
// forever with no error).
def optional_file_ch(paramValue, sentinel) {
    return paramValue ?
        channel.fromPath(paramValue, type: 'file', checkIfExists: true).first() :
        channel.value(file(sentinel))
}

// Collect the files of a bowtie2 (.bt2/.bt2l) or hisat2 (.ht2/.ht2l) index
// from its basename into a value channel. The files are passed to the
// aligners as `path` inputs, so Nextflow stages them into each task and makes
// them visible inside the container wherever they live on disk.
def index_files_ch(basename, ext) {
    def idx = files("${basename}.*.${ext}*")
    assert idx.size() > 0 : "No ${ext} index files found for basename '${basename}' (expected ${basename}.1.${ext} etc.)"
    return channel.value(idx)
}

// Memory of the SLURM job this process runs inside (null when not in a job).
// SLURM_MEM_PER_NODE is set by --mem, SLURM_MEM_PER_CPU by --mem-per-cpu;
// both are in MB. Mirrors the calculation used by -profile standard in
// nextflow.config.
def slurm_job_memory() {
    def env = System.getenv()
    def cpus = (env.SLURM_CPUS_PER_TASK ?: env.SLURM_CPUS_ON_NODE ?: '1') as long
    if (env.SLURM_MEM_PER_NODE) { return MemoryUnit.of("${env.SLURM_MEM_PER_NODE} MB") }
    if (env.SLURM_MEM_PER_CPU)  { return MemoryUnit.of("${(env.SLURM_MEM_PER_CPU as long) * cpus} MB") }
    return null
}

// Which reference a BAM was aligned to: BOWTIE2 names off-target
// alignments "<sample>_ecoli" / "<sample>_spikein" (simpleName keeps
// everything before the first '.', i.e. before .MAPQ*.NoDups.bam).
def reference_of(bam) {
    if (bam.simpleName.endsWith('_ecoli'))   { return 'ecoli' }
    if (bam.simpleName.endsWith('_spikein')) { return 'spikein' }
    return 'genome'
}

// Write <reads dir>/QC/<sample>_commands.txt (next to the other QC files): every
// command the pipeline ran for that sample, in the order the tasks started.
// Reads the run's trace file (one row per task, with its work directory) and
// each task's .command.sh, the exact script Nextflow executed. Called at the
// end of the run, also when it fails, so the file shows how far the sample
// got and the command that failed. On -resume, cached tasks are included with
// their original start times. Tasks are matched to samples by the tag prefix
// "<sample dir name> | ..." used by every process in modules.nf.
def write_command_logs() {
    def traceFile = file("${(params.outdir as String) ?: (params.target as String) ?: '.'}/pipeline_info/nf-trace_${params.run_timestamp}.txt")
    if (!traceFile.exists()) {
        log.warn "Trace file not found (${traceFile}) - per-sample command logs not written"
        return
    }
    def lines = traceFile.readLines()
    def header = lines[0].tokenize('\t')
    def rows = lines.drop(1).findAll { String l -> l.trim() }.collect { String l ->
        def cols = l.split('\t', -1) as List
        def row = [:]
        header.eachWithIndex { String h, int i -> row[h] = (i < cols.size()) ? cols[i] : '' }
        row
    }

    // Written next to the other QC files, in each sample's reads directory.
    def sampleDirs = [file(params.target as String)]
    if (params.control) {
        sampleDirs << file(params.control as String)
    }

    sampleDirs.each { sampleDir ->
        def sampleName = sampleDir.name
        def tasks = rows.findAll { Map r -> (r.tag as String).tokenize('|')[0]?.trim() == sampleName }
        // Execution order: start time (sortable text), then task id; tasks that never started go last
        tasks = tasks.sort { Map r ->
            def start = (r.start as String) in ['', '-'] ? '9999' : (r.start as String)
            "${start} ${(r.task_id as String).padLeft(8, '0')}"
        }
        def out = new StringBuilder()
        out << "# Commands run by ${workflow.manifest.name} ${workflow.manifest.version} for sample: ${sampleName}\n"
        out << "# Run:       ${workflow.runName} (session ${workflow.sessionId})\n"
        out << "# Launched:  ${workflow.commandLine}\n"
        out << "# Nextflow:  ${nextflow.version}   Profile: ${workflow.profile}\n"
        if (workflow.containerEngine) {
            out << "# Container: ${workflow.containerEngine} ${workflow.container}\n"
        }
        out << "# Finished:  ${workflow.complete}   Status: ${workflow.success ? 'success' : 'FAILED'}\n"
        out << "# Each block below is one task, in the order it started. CACHED = reused from an\n"
        out << "# earlier run (-resume). File names are relative to that task's work directory.\n"
        tasks.eachWithIndex { Map r, int i ->
            def script = file("${r.workdir}/.command.sh")
            def body = script.exists() ?
                script.readLines().findAll { String l -> !l.startsWith('#!') }.join('\n').trim() :
                '(script not found - work directory removed?)'
            out << "\n# " + ('-' * 76) + "\n"
            out << "# [${i + 1}] ${r.name}   status=${r.status} exit=${r.exit} attempt=${r.attempt}\n"
            out << "# started ${r.start}   work dir: ${r.workdir}\n"
            out << "# " + ('-' * 76) + "\n"
            out << body << "\n"
        }
        def qcDir = sampleDir.resolve('QC')
        qcDir.mkdirs()
        def dest = qcDir.resolve("${sampleName}_commands.txt")
        dest.text = out.toString()
        log.info "Commands for ${sampleName} (${tasks.size()} tasks) written to ${dest}"
    }
}

workflow {

    display_start()

    assert params.target       : "params.target is required (directory of target fastq.gz reads)"
    assert params.genome_index : "params.genome_index is required (bowtie2/hisat2 index basename)"

    // Running locally (-profile standard) inside a SLURM job: every task
    // shares that job's memory, which the local executor cannot enforce per
    // task. Stop now with a clear message rather than have the aligner killed
    // (exit 137) part-way through.
    def jobMem = slurm_job_memory()
    if (jobMem && !workflow.profile.tokenize(',').contains('slurm')) {
        def minMem = params.min_memory as MemoryUnit
        assert jobMem >= minMem :
            "This SLURM job has ${jobMem} of memory, but running the pipeline inside one job needs at least ${minMem} (the mm10 bowtie2 index alone needs ~4 GB). Request more, e.g. --mem-per-cpu=2G or --mem=64G, or lower --min_memory"
        log.info "Running inside SLURM job ${System.getenv('SLURM_JOB_ID')}: ${jobMem} memory available"
    }

    // Samples are identified by their reads directory, which travels
    // unchanged as the first element of every tuple. Branching compares that
    // directory name exactly, so it works whatever the directory is called
    // (dots, or a control name that contains the target name).
    def merge_mode = params.merge && params.reads_type != 'paired'
    if (params.merge && !merge_mode) {
        log.warn "--merge only applies to single-end reads (--reads_type single) and is ignored for paired-end data"
    }

    def targetName  = file(params.target).name
    def controlName = params.control ? file(params.control).name : null
    if (controlName) {
        assert targetName != controlName : "Target and control directories must have different names (both are '${targetName}')"
        // Output files are named with the directory name up to its first dot,
        // so those prefixes must differ too or the two samples' files collide.
        assert file(params.target).simpleName != file(params.control).simpleName :
            "Target '${targetName}' and control '${controlName}' share the prefix '${file(params.target).simpleName}' (text before the first '.'), so their output files would collide. Rename one directory"
    }

    target_reads  = get_reads_ch( params.target )
    control_reads = get_reads_ch( params.control )
    raw_reads     = target_reads.concat( control_reads )

    // One MD5SUMCHECK task per individual fastq.gz file
    raw_reads
        .flatMap { sample, reads -> (reads instanceof List ? reads : [reads]).collect { r -> tuple(sample, r) } }
        .set { individual_reads }
    MD5SUMCHECK( individual_reads )

    // Reference indexes: (files, basename) for each aligner call
    def genome_ext  = (params.assay == 'rnaseq') ? 'ht2' : 'bt2'
    genome_files    = index_files_ch( params.genome_index, genome_ext )
    genome_name     = file(params.genome_index).name

    blacklist_regions = optional_file_ch( params.blacklist,      'NO_BLACKLIST_FILE' )
    multiqc_config     = optional_file_ch( params.multiqc_config, 'NO_MULTIQC_CONFIG_FILE' )

    def hg38_egs = [50: '2701495761', 75: '2747877777', 100: '2805636331', 150: '2862010578', 200: '2887553303']
    def mm10_egs = [50: '2308125349', 75: '2407883318', 100: '2467481108', 150: '2494787188', 200: '2520869189']
    def rn6_egs  = [50: '2375372135', 75: '2440746491', 100: '2480029900', 150: '2477334634', 200: '2478552171']
    def egs_by_assembly = [mm10: mm10_egs, hg38: hg38_egs, rn6: rn6_egs]
    assert egs_by_assembly.containsKey(params.assembly) : "Assembly '${params.assembly}' not available to estimate effective_genome_size (choose from: ${egs_by_assembly.keySet()})"
    assert egs_by_assembly[params.assembly].containsKey(params.length as int) : "No effective_genome_size entry for read length ${params.length} (choose from: ${egs_by_assembly[params.assembly].keySet()})"
    effective_genome_size = channel.value( egs_by_assembly[params.assembly][params.length as int] )

    // File integrity check already runs on individual_reads above.

    // FASTQC on original reads
    FASTQC_RAW( raw_reads )

    // Trim if facility provided as raw
    if (!params.pretrimmed) {
        CUTADAPT( raw_reads )
        FASTQC_TRIMMED( CUTADAPT.out.trimmed_reads )
        reads = CUTADAPT.out.trimmed_reads
    } else {
        reads = raw_reads
    }

    // Align reads
    if (params.assay == 'rnaseq') {

        ALIGN_RNA_MAIN( reads, genome_files, genome_name )
        initial_bams = ALIGN_RNA_MAIN.out.mapped_reads

    } else {

        ALIGN_MAIN( reads, genome_files, genome_name, channel.value('target') )
        initial_bams = ALIGN_MAIN.out.mapped_reads

        if (params.ecoli_index) {
            ALIGN_ECOLI( reads, index_files_ch(params.ecoli_index, 'bt2'), file(params.ecoli_index).name, channel.value('ecoli') )
            initial_bams = initial_bams.concat( ALIGN_ECOLI.out.mapped_reads )
        }
        if (params.spikein && params.spikein_index) {
            ALIGN_SPIKEIN( reads, index_files_ch(params.spikein_index, 'bt2'), file(params.spikein_index).name, channel.value('spikein') )
            initial_bams = initial_bams.concat( ALIGN_SPIKEIN.out.mapped_reads )
        }
    }

    // Alignments arrive coordinate-sorted from the aligners. Pre-filter QC
    // (flagstat + duplication rate) runs alongside the filtering steps.
    PREFILTER_QC( initial_bams )

    // Filter out poor alignments (indexed in the same step)
    QFILTER_BAM( initial_bams )

    // Remove duplicates (indexed in the same step)
    DEDUPLICATE_BAM( QFILTER_BAM.out.bam_qfiltered )

    // Final alignments: both the MAPQ-filtered and the de-duplicated BAM are kept
    indexed_bams = QFILTER_BAM.out.bam_qfiltered.mix( DEDUPLICATE_BAM.out.bam_deduplicated )

    // --merge (single-end): up to here each fastq.gz was processed as its own
    // library. Merge each sample's per-file BAMs - separately for each reference
    // (genome / ecoli / spikein) and for the MAPQ-filtered and de-duplicated
    // versions - into the sample's BAM, named as without --merge. Everything
    // below runs on the merged BAMs; the per-file BAMs also get their own bigwigs.
    if (merge_mode) {
        // number of files per sample, so each merge starts as soon as its last file is ready
        def fileCounts = [params.target, params.control]
            .findAll { d -> d instanceof CharSequence && d }
            .collectEntries { d -> [(file(d as String).name): files("${d}/*.fastq.gz").size()] }
        indexed_bams
            .map { sample, bam, bai ->
                def ref = reference_of(bam)
                def mergedName = sample.simpleName + (ref == 'genome' ? '' : "_${ref}") + bam.name.substring(bam.name.indexOf('.'))
                tuple(groupKey("${sample.name}/${mergedName}", fileCounts[sample.name] as int), sample, mergedName, bam, bai)
            }
            .groupTuple()
            .map { _key, samples, names, bams, bais -> tuple(samples[0], names[0], bams, bais) }
            .set { bams_to_merge }
        MERGE_BAMS( bams_to_merge )
        final_bams     = MERGE_BAMS.out.bam_merged
        replicate_bams = indexed_bams
    } else {
        final_bams     = indexed_bams
        replicate_bams = channel.empty()
    }

    if (params.assay == 'rnaseq') {

        final_bams
            .branch { sample, bam, _bai ->
                def isDedup = bam.name.contains('.NoDups.')
                target:  sample.name == targetName  && !isDedup
                control: controlName && sample.name == controlName && !isDedup
            }.set { bams }

        BIGWIGS_TARGET_FORWARD( bams.target, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('forward') )
        BIGWIGS_TARGET_REVERSE( bams.target, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('reverse') )
        // Both strands combined (all reads)
        BIGWIGS_TARGET_UNSTRANDED( bams.target, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('unstranded') )
        if (params.control) {
            BIGWIGS_CONTROL_FORWARD( bams.control, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('forward') )
            BIGWIGS_CONTROL_REVERSE( bams.control, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('reverse') )
            BIGWIGS_CONTROL_UNSTRANDED( bams.control, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('unstranded') )
        }
        if (merge_mode) {
            // per-file bigwigs (same strand set as the merged BAM)
            replicate_bams
                .filter { _sample, bam, _bai -> reference_of(bam) == 'genome' && !bam.name.contains('.NoDups.') }
                .set { rep_bams }
            BIGWIGS_REPLICATES_FORWARD( rep_bams, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('forward') )
            BIGWIGS_REPLICATES_REVERSE( rep_bams, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('reverse') )
            BIGWIGS_REPLICATES_UNSTRANDED( rep_bams, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('unstranded') )
        }

        MULTIQC( params.control ? bams.target.concat(bams.control) : bams.target, multiqc_config )

    } else {

        final_bams
            .branch { sample, bam, _bai ->
                def isDedup   = bam.name.contains('.NoDups.')
                def isTarget  = sample.name == targetName
                def isControl = controlName && sample.name == controlName
                def ref       = reference_of(bam)
                target:          isTarget  && ref == 'genome'  && isDedup
                target_ecoli:    isTarget  && ref == 'ecoli'   && isDedup
                target_spikein:  isTarget  && ref == 'spikein' && isDedup
                control:         isControl && ref == 'genome'  && isDedup
                control_ecoli:   isControl && ref == 'ecoli'   && isDedup
                control_spikein: isControl && ref == 'spikein' && isDedup
            }.set { bams }

        if (params.ecoli_index) {
            ECOLI_SCALE_FACTOR_TARGET( bams.target, bams.target_ecoli )
            if (params.control) {
                ECOLI_SCALE_FACTOR_CONTROL( bams.control, bams.control_ecoli )
            }
        }

        // Bigwigs: unscaled and ecoli-normalized
        BIGWIGS_NO_SCALE_FACTOR_TARGET( bams.target, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('na') )
        if (merge_mode) {
            // per-file bigwigs (unscaled, from each file's de-duplicated BAM)
            BIGWIGS_REPLICATES(
                replicate_bams.filter { _sample, bam, _bai -> reference_of(bam) == 'genome' && bam.name.contains('.NoDups.') },
                blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('na') )
        }
        if (params.ecoli_index) {
            BIGWIGS_ECOLI_SCALE_FACTOR_TARGET( bams.target, blacklist_regions, effective_genome_size, ECOLI_SCALE_FACTOR_TARGET.out.scale_factor, channel.value('na') )
        }
        if (params.control) {
            BIGWIGS_NO_SCALE_FACTOR_CONTROL( bams.control, blacklist_regions, effective_genome_size, channel.value(tuple('None', 1)), channel.value('na') )
            if (params.ecoli_index) {
                BIGWIGS_ECOLI_SCALE_FACTOR_CONTROL( bams.control, blacklist_regions, effective_genome_size, ECOLI_SCALE_FACTOR_CONTROL.out.scale_factor, channel.value('na') )
            }
            BIGWIG_BAMCOMPARE( bams.target, bams.control, blacklist_regions, effective_genome_size, channel.value( file("${params.genome_index}.fa", checkIfExists: true) ) )
        }

        // Bigwigs: spike-in normalized
        if (params.spikein && params.spikein_index) {
            SPIKEIN_SCALE_FACTOR_TARGET( bams.target, bams.target_spikein )
            BIGWIGS_SPIKEIN_SCALE_FACTOR_TARGET( bams.target, blacklist_regions, effective_genome_size, SPIKEIN_SCALE_FACTOR_TARGET.out.scale_factor, channel.value('na') )
            if (params.control) {
                SPIKEIN_SCALE_FACTOR_CONTROL( bams.control, bams.control_spikein )
                BIGWIGS_SPIKEIN_SCALE_FACTOR_CONTROL( bams.control, blacklist_regions, effective_genome_size, SPIKEIN_SCALE_FACTOR_CONTROL.out.scale_factor, channel.value('na') )
            }
        }

        // Peak calling
        def targetControlBams = params.control ? bams.target.concat(bams.control) : bams.target

        BEDGRAPH_COVERAGE( targetControlBams, blacklist_regions, effective_genome_size )
        BEDGRAPH_COVERAGE.out.bedgraph
            .branch { sample, _bedgraph ->
                target:  sample.name == targetName
                control: controlName && sample.name == controlName
            }.set { bedgraphs }

        MACS_PEAKS( targetControlBams, effective_genome_size )
        SEACR_PEAKS( params.control ? bedgraphs.target.concat(bedgraphs.control) : bedgraphs.target )
        GOPEAKS_PEAKS( targetControlBams )

        // Genrich needs BAMs sorted by read name: sort each BAM once and share
        // it between the per-sample and target-vs-control calls
        NAME_SORT_BAM( targetControlBams )
        NAME_SORT_BAM.out.bam_namesorted
            .branch { sample, _bam ->
                target:  sample.name == targetName
                control: controlName && sample.name == controlName
            }.set { namesorted_bams }
        GENRICH_PEAKS( NAME_SORT_BAM.out.bam_namesorted, blacklist_regions )

        if (params.control) {
            MACS_PEAKS_WITH_CONTROL( bams.target, bams.control, effective_genome_size )
            SEACR_PEAKS_WITH_CONTROL( bedgraphs.target, bedgraphs.control )
            GOPEAKS_PEAKS_WITH_CONTROL( bams.target, bams.control )
            GENRICH_PEAKS_WITH_CONTROL( namesorted_bams.target, namesorted_bams.control, blacklist_regions )
        }

        // FRiP: gather every peak file each sample produced (peak callers that
        // failed and were skipped simply contribute nothing) and score them all
        // against that sample's final BAM in one task per sample
        peak_files = MACS_PEAKS.out.macs_narrowpeak
            .mix( MACS_PEAKS.out.macs_broadpeak, MACS_PEAKS.out.macs_gappedpeak,
                  SEACR_PEAKS.out.seacr_peaks, GOPEAKS_PEAKS.out.gopeaks_peaks,
                  GENRICH_PEAKS.out.genrich_peaks )
        if (params.control) {
            peak_files = peak_files.mix(
                MACS_PEAKS_WITH_CONTROL.out.macs_narrowpeak, MACS_PEAKS_WITH_CONTROL.out.macs_broadpeak,
                MACS_PEAKS_WITH_CONTROL.out.macs_gappedpeak, SEACR_PEAKS_WITH_CONTROL.out.seacr_peaks,
                GOPEAKS_PEAKS_WITH_CONTROL.out.gopeaks_peaks, GENRICH_PEAKS_WITH_CONTROL.out.genrich_peaks )
        }
        targetControlBams
            .map { sample, bam, bai -> tuple(sample.name, sample, bam, bai) }
            .join( peak_files.groupTuple() )
            .map { _name, sample, bam, bai, peaks -> tuple(sample, bam, bai, peaks) }
            .set { frip_input }
        FRIP( frip_input )

        MULTIQC( targetControlBams, multiqc_config )
    }

    // Assignment form inside the entry workflow: accepted by both the
    // classic parser and the strict syntax parser (Nextflow 25.10+ default),
    // unlike a top-level `workflow.onComplete { }` block.
    workflow.onComplete = {
        try {
            write_command_logs()
        } catch (Exception e) {
            log.warn "Could not write per-sample command logs: ${e.message}"
        }
        log.info "Pipeline ${workflow.success ? 'completed successfully' : 'failed'} at ${workflow.complete} (duration: ${workflow.duration})"
    }
}

