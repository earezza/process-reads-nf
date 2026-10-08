/*
 * ===================================================
 *  Process Reads Pipeline - process definitions
 * ===================================================
 * Resource requests (cpus/memory/time) are intentionally NOT hardcoded here -
 * they are assigned centrally per label in nextflow.config so the same
 * modules run efficiently on a laptop (-profile standard) or a SLURM cluster
 * (-profile slurm) without editing this file.
 */

// Groovy helper: Nextflow presents a `path` input as a single Path when only
// one file is staged, and as a List<Path> when several are. Normalize to a
// List so processes behave the same for paired-end (2 files) and single-end
// (1 file) reads.
def asList(x) { return (x instanceof List) ? x : [x] }

// Each sample's original reads directory (the --target or --control path), so
// every output is published next to the raw reads. A process's `sample` input is
// the staged copy of that directory, which on its own carries only the folder
// name - used directly in publishDir it would resolve against the launch
// directory and create a second folder there.
def readsDir(sample) {
    def match = [params.target, params.control]
        .findAll { d -> d instanceof CharSequence && d }
        .collect { d -> file(d as String) }
        .find { d -> d.name == sample.name }
    return match ? match.toString() : sample.name
}

// --merge with single-end reads: every fastq.gz in a sample's directory is a
// separate library, processed on its own up to its final BAM; the per-file BAMs
// are then merged into the sample's BAM (MERGE_BAMS). Per-file outputs are
// named after the file (rep1.fastq.gz -> rep1); otherwise after the sample
// directory, as without --merge.
def isMerge() {
    return params.merge && params.reads_type != 'paired'
}
def readsId(f) {
    return f.name.replaceAll('(_trimmed)?\\.(fastq|fq)(\\.gz)?$', '')
}
def unitName(sample, reads) {
    return isMerge() ? readsId(asList(reads)[0]) : sample.simpleName
}
def unitTag(sample, reads) {
    return isMerge() ? "${sample.name} | ${readsId(asList(reads)[0])}" : "${sample.name}"
}

// Read-group fields. In merge mode each file gets its own read-group ID and
// library (LB) and shares the sample name (SM), so they stay distinguishable
// in the merged BAM.
def readGroupArgs(sample, unit) {
    def rg_id = isMerge() ? (params.rg_id ? "${params.rg_id}.${unit}" : unit) : (params.rg_id ?: sample.name)
    def rg_sm = params.rg_sm ?: (isMerge() ? sample.name : rg_id)
    def rg_lb = params.rg_lb ?: (isMerge() ? unit : params.reads_type)
    return "--rg-id ${rg_id} --rg SM:${rg_sm} --rg LB:${rg_lb} --rg PU:${params.rg_pu} --rg PL:${params.rg_pl}"
}

// Shell snippet: run an aligner piped into samtools sort and exit with the status
// of whichever side actually failed. The bowtie2/hisat2 wrapper scripts exit 1
// even when the aligner binary was killed (e.g. 137 = out of memory) and only
// write the real code into the log, and with `pipefail` alone the task would
// report samtools' secondary "fail to read the header" error instead. Passing
// the real code through lets errorStrategy retry OOM kills with more memory;
// printing the log puts the aligner's own message in Nextflow's error report.
// (Lines are indented to match the process script blocks they are inserted
// into, so Nextflow's indent stripping keeps `#!/bin/bash` in column 1.)
def alignPipe(align_cmd, log_file, bam_file, cpus) {
    // Alignments are coordinate-sorted as they stream out of the aligner (no
    // separate sort step or unsorted BAM on disk). samtools sort mostly waits
    // while the aligner runs and uses its threads for the final merge.
    def sort_threads = Math.max(1, (cpus as int).intdiv(2))
    def sort_order   = (params.sort_bam_order == 'queryname') ? '-n ' : ''
    return [
        "set +e",
        "${align_cmd} 2> ${log_file} | samtools sort -@ ${sort_threads} ${sort_order}-T ${bam_file}.tmp -o ${bam_file} -",
        'status=("${PIPESTATUS[@]}")',
        "set -e",
        'if [ "${status[0]}" -ne 0 ]; then',
        '    echo "Aligner failed (exit ${status[0]}); log follows:" >&2',
        "    cat ${log_file} >&2",
        "    code=\$(sed -n 's/.*exited with value \\([0-9][0-9]*\\).*/\\1/p' ${log_file} | tail -n 1)",
        '    exit "${code:-${status[0]}}"',
        "fi",
        'exit "${status[1]}"'
    ].join('\n        ')
}

// Java heap for Picard, sized from the task's memory request (85%, leaving room
// for the JVM itself). Without it Java picks a default - typically a quarter of
// the memory it can see, which is either too small for large BAMs or, under the
// local executor, more than the task asked for.
def javaHeap(mem) {
    return "-Xmx${(mem.toMega() * 85).intdiv(100)}m"
}


process MD5SUMCHECK {
    label 'md5sum_check'
    tag "${sample.name} | ${read}"
    publishDir { "${readsDir(sample)}" }, mode: params.publish_mode, overwrite: false, pattern: "*.md5"

    input:
        tuple path(sample), path(read)

    output:
        path "${read}.md5"

    script:
        // The read is staged directly in the work dir; a pre-existing hash
        // (if any) lives alongside the original file inside the staged
        // `sample` directory. Checking with a shell test avoids declaring a
        // Nextflow `path` input for a file that may not exist.
        """
        #!/bin/sh
        set -e
        if [ -f "${sample}/${read}.md5" ]; then
            echo "Found existing md5 checksum for ${read}, verifying..."
            cp "${sample}/${read}.md5" "${read}.md5"
        else
            echo "Generating md5 checksum hash for ${read}..."
            md5sum ${read} > ${read}.md5
        fi
        md5sum -c ${read}.md5
        """
}

process FASTQC {
    label 'fastqc_check'
    tag "${unitTag(sample, reads)}"
    cpus { Math.min(asList(reads).size(), 2) }    // FastQC runs one thread per file
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*fastqc.{zip,html}"

    input:
        tuple path(sample), path(reads)

    output:
        path "*fastqc.{zip,html}"

    script:
        """
        #!/bin/bash
        set -euo pipefail
        fastqc -t ${task.cpus} ${reads}
        """
}

process MULTIQC {
    label 'multiqc_report'
    tag "${sample.name}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*.html"

    input:
        tuple path(sample), path(bam), path(bai)
        path multiqc_config

    output:
        path "${sample.simpleName}_multiqc_report.html"

    script:
        config_file = (multiqc_config.name == 'NO_MULTIQC_CONFIG_FILE') ? '' : "-c ${multiqc_config}"
        """
        #!/bin/bash
        set -euo pipefail
        multiqc ${sample} ${config_file} -v --force --filename ${sample.simpleName}_multiqc_report.html
        """
}

// Adapter / quality trimming in one cutadapt pass per file (RNA-seq poly-A
// trimming included - cutadapt applies it after adapter trimming). Output is
// gzip level 1 (-Z): written in parallel by cutadapt's cores and several times
// smaller than plain FASTQ, which cuts disk I/O for this and the next steps.
process CUTADAPT {
    label 'cutadapt'
    tag "${unitTag(sample, reads)}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "cutadapt_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        path "cutadapt_${unit}.log"
        tuple path(sample), path("*_trimmed.fastq.gz"), emit: trimmed_reads

    script:
        def r = asList(reads)
        unit = unitName(sample, reads)
        polya = (params.assay == 'rnaseq' && params.polyAtrim) ? '--poly-a' : ''
        if (params.reads_type == 'paired') {
            adapters = params.trim_adapters ? '-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT' : ''
            trim     = params.qctrim ? '-m 20 -q 20,20 -u 11 -U 11' : ''
            outputs  = "-o ${r[0].name.replace('.fastq.gz', '_trimmed.fastq.gz')} -p ${r[1].name.replace('.fastq.gz', '_trimmed.fastq.gz')}"
            """
            #!/bin/bash
            set -euo pipefail
            cutadapt --cores=${task.cpus} -Z ${trim} ${adapters} ${polya} ${outputs} ${r[0]} ${r[1]} &>> cutadapt_${unit}.log
            """
        } else {
            // Single-end: one file (in merge mode each file is its own task)
            adapters = params.trim_adapters ? '-a CTGTCTCTTATACACATCT' : ''
            trim     = params.qctrim ? '-m 20 -q 20 -u 11' : ''
            """
            #!/bin/bash
            set -euo pipefail
            for f in ${r.join(' ')}; do
                base=\$(basename "\$f" .fastq.gz)
                cutadapt --cores=${task.cpus} -Z ${trim} ${adapters} ${polya} -o "\${base}_trimmed.fastq.gz" "\$f" &>> cutadapt_${unit}.log
            done
            """
        }
}


process BOWTIE2 {
    label 'bowtie2'
    tag { "${unitTag(sample, reads)} | ${align_to}" }
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "bowtie2_*.log"

    input:
        tuple path(sample), path(reads)
        // Index files are staged as inputs (not passed as a path string) so
        // Nextflow makes them visible inside the container on any filesystem.
        path(index_files, stageAs: 'index/*')
        val(index_name)
        val(align_to)

    output:
        path "bowtie2_${unit}_${align_to}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads

    script:
        def r = asList(reads)
        unit = unitName(sample, reads)
        reads_name = (align_to in ['target', 'control']) ? "${unit}" : "${unit}_${align_to}"

        if (params.reads_type == 'paired') {
            align = (align_to == 'spikein') ?
                "-1 ${r[0]} -2 ${r[1]} --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -k 1" :
                (align_to == 'ecoli') ?
                "-1 ${r[0]} -2 ${r[1]} --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700" :
                "-1 ${r[0]} -2 ${r[1]} --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
        } else {
            reads_csv = r.join(',')
            align = (align_to == 'spikein') ?
                "-U ${reads_csv} --end-to-end --very-sensitive --phred33 -k 1" :
                (align_to == 'ecoli') ?
                "-U ${reads_csv} --end-to-end --very-sensitive --no-unal --phred33" :
                "-U ${reads_csv} --local --very-sensitive-local --no-unal --phred33"
        }

        align_cmd = "bowtie2 -p ${task.cpus} ${align} -x index/${index_name} ${readGroupArgs(sample, unit)}"
        """
        #!/bin/bash
        set -euo pipefail
        ${alignPipe(align_cmd, "bowtie2_${unit}_${align_to}.log", "${reads_name}.bam", task.cpus)}
        """
}

process HISAT2 {
    label 'hisat2'
    tag "${unitTag(sample, reads)}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "hisat2_*.log"

    input:
        tuple path(sample), path(reads)
        path(index_files, stageAs: 'index/*')
        val(index_name)

    output:
        path "hisat2_${unit}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads

    script:
        def r = asList(reads)
        unit = unitName(sample, reads)
        reads_name = "${unit}"
        // The aligner binary (hisat2-align-s/-l) is run directly rather than via
        // the `hisat2` launcher script. Given .gz files, the launcher decompresses
        // them through named pipes (or, with --no-named-pipes, full temporary
        // files) called /tmp/<process ID>.*; inside the container every task has
        // the same small process ID while /tmp is the node's, so alignments
        // running side by side collide ("mkfifo(/tmp/42.unp) failed"). Passing
        // decompressed streams instead fails the launcher's own regular-file
        // check ("Read file '/dev/fd/63' doesn't exist"). The binary has neither
        // problem: it reads the streams below, decompressed by separate gzip
        // processes alongside the alignment threads. `--wrapper basic-0` is what
        // the launcher itself passes; options and log output are unchanged.
        align = (params.reads_type == 'paired') ?
            "-1 <(gzip -cdf ${r[0]}) -2 <(gzip -cdf ${r[1]}) --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700" :
            "-U <(gzip -cdf ${r.join(' ')}) --no-unal --phred33"

        align_cmd = "\"\$hisat2_bin\" --wrapper basic-0 -p ${task.cpus} ${align} -x index/${index_name} ${readGroupArgs(sample, unit)}"
        """
        #!/bin/bash
        set -euo pipefail
        # Locate the aligner binary next to the launcher (as the launcher does),
        # using the large-index build when only a .ht2l index is present.
        hisat2_dir=\$(dirname "\$(readlink -f "\$(command -v hisat2)")")
        hisat2_bin="\$hisat2_dir/hisat2-align-s"
        if [ ! -e index/${index_name}.1.ht2 ] && [ -e index/${index_name}.1.ht2l ]; then
            hisat2_bin="\$hisat2_dir/hisat2-align-l"
        fi
        ${alignPipe(align_cmd, "hisat2_${unit}.log", "${reads_name}.bam", task.cpus)}
        """
}

// Pre-filter QC on each sorted alignment: read counts and duplication rate before
// MAPQ filtering. Runs alongside the filtering steps, not in front of them.
process PREFILTER_QC {
    label 'prefilter_qc'
    tag "${sample.name} | ${bam.name}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*.log"

    input:
        tuple path(sample), path(bam)

    output:
        path "picard-dupStats_${bam.simpleName}.log"
        path "samtools-flagstat-prefiltered_${bam.simpleName}.log"

    script:
        """
        #!/bin/bash
        set -euo pipefail
        export JAVA_TOOL_OPTIONS="${javaHeap(task.memory)}"
        samtools flagstat -@ ${task.cpus} ${bam} &>> samtools-flagstat-prefiltered_${bam.simpleName}.log
        # Duplicate stats only (duplicates are removed after MAPQ filtering, in DEDUPLICATE_BAM)
        picard MarkDuplicates -I ${bam} -O /dev/null -METRICS_FILE picard-dupStats_${bam.simpleName}.log
        """
}

// MAPQ / pairing filter. The input is already coordinate-sorted and filtering
// keeps that order, so no re-sort is needed; the result is indexed here.
process QFILTER_BAM {
    label 'qfilter_bam'
    tag "${sample.name} | ${bam.name}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*.log"
    publishDir { "${readsDir(sample)}/bams/" }, mode: params.publish_mode, overwrite: false, pattern: "*.{bam,bai}"

    input:
        tuple path(sample), path(bam)

    output:
        path "samtools-flagstat-postqfiltered_${bam.simpleName}.log"
        tuple path(sample), path("${out_bam}"), path("${out_bam}.bai"), emit: bam_qfiltered

    script:
        sam_flags = params.reads_type == 'paired' ? '-f 2' : '-F 4'
        out_bam = "${bam.simpleName}.MAPQ${params.mapq}.bam"
        """
        #!/bin/bash
        set -euo pipefail
        samtools view -@ ${task.cpus} -b ${sam_flags} -q ${params.mapq} -o ${out_bam} ${bam}
        samtools index -@ ${task.cpus} ${out_bam}
        samtools flagstat -@ ${task.cpus} ${out_bam} &>> samtools-flagstat-postqfiltered_${bam.simpleName}.log
        """
}

process DEDUPLICATE_BAM {
    label 'deduplicate_bam'
    tag "${sample.name} | ${bam.name}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*.log"
    publishDir { "${readsDir(sample)}/bams/" }, mode: params.publish_mode, overwrite: false, pattern: "*.{bam,bai}"

    input:
        tuple path(sample), path(bam), path(bam_index)

    output:
        path "picard-deduplicate_${bam.simpleName}.log"
        path "samtools-flagstat_${bam.simpleName}.log"
        tuple path(sample), path("${bam.baseName}.NoDups.bam"), path("${bam.baseName}.NoDups.bam.bai"), emit: bam_deduplicated

    script:
        """
        #!/bin/bash
        set -euo pipefail
        export JAVA_TOOL_OPTIONS="${javaHeap(task.memory)}"
        picard MarkDuplicates -I ${bam} -O ${bam.baseName}.NoDups.bam -REMOVE_DUPLICATES true -METRICS_FILE picard-deduplicate_${bam.simpleName}.log
        samtools index -@ ${task.cpus} ${bam.baseName}.NoDups.bam
        samtools flagstat -@ ${task.cpus} ${bam.baseName}.NoDups.bam &>> samtools-flagstat_${bam.simpleName}.log
        """
}

// Merge mode: combine one sample's per-file final BAMs (each already
// MAPQ-filtered and de-duplicated as its own library) into the sample's BAM,
// named as it would be without --merge so every later step is unchanged.
process MERGE_BAMS {
    label 'merge_bams'
    tag "${sample.name} | ${merged_name}"
    publishDir { "${readsDir(sample)}/bams/" }, mode: params.publish_mode, overwrite: false, pattern: "*.{bam,bai}"

    input:
        tuple path(sample), val(merged_name), path(bams, stageAs: 'in/*'), path(bam_index, stageAs: 'in/*')

    output:
        tuple path(sample), path("${merged_name}"), path("${merged_name}.bai"), emit: bam_merged

    script:
        def input_args = asList(bams).collect { b -> "-I ${b}" }.join(' ')
        """
        #!/bin/bash
        set -euo pipefail
        export JAVA_TOOL_OPTIONS="${javaHeap(task.memory)}"
        picard MergeSamFiles ${input_args} -O ${merged_name} -SORT_ORDER coordinate -USE_THREADING true
        samtools index -@ ${task.cpus} ${merged_name}
        """
}

process GET_SCALE_FACTOR {
    label 'scale_factor'
    tag "${sample_target.name} | ${bam_norm.simpleName}"
    publishDir { "${readsDir(sample_target)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*.txt"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        // stageAs avoids a work-dir name collision: sample_norm is the same
        // physical reads directory as sample_target (e.g. the target's own
        // ecoli/spike-in alignment), so it must be staged under a distinct name.
        tuple path(sample_norm, stageAs: 'sample_norm'), path(bam_norm), path(bam_index_norm)

    output:
        tuple val("${bam_norm.simpleName}"), stdout, emit: scale_factor
        path "scale_factor_${bam_target.simpleName}_${bam_norm.simpleName}.txt"

    script:
        """
        #!/bin/bash
        set -euo pipefail
        # Read counts from the BAM indexes (idxstats: mapped + unmapped per reference),
        # the same totals as samtools view -c without reading the whole BAMs
        count_target=\$(samtools idxstats ${bam_target} | awk '{s += \$3 + \$4} END {print s}')
        count_norm=\$(samtools idxstats ${bam_norm} | awk '{s += \$3 + \$4} END {print s}')
        echo "${bam_target}: \$count_target" > scale_factor_${bam_target.simpleName}_${bam_norm.simpleName}.txt
        echo "${bam_norm}: \$count_norm" >> scale_factor_${bam_target.simpleName}_${bam_norm.simpleName}.txt
        awk -v t="\$count_target" -v s="\$count_norm" 'BEGIN { print t / s }' | tee -a scale_factor_${bam_target.simpleName}_${bam_norm.simpleName}.txt
        """
}

process BIGWIG_COVERAGE {
    label 'bigwig'
    tag "${sample.name} | ${bam.simpleName}"
    publishDir { "${readsDir(sample)}/bigwigs/" }, mode: params.publish_mode, overwrite: false, pattern: "*.bw"

    input:
        tuple path(sample), path(bam), path(bam_index)
        path(blacklist_file)
        val(effective_genome_size)
        tuple val(scale_source), val(scale_factor)
        val(direction)

    output:
        path "${bw_name}", emit: bigwig

    script:
        blacklist = (blacklist_file.name == 'NO_BLACKLIST_FILE') ? '' : "--blackListFileName ${blacklist_file}"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' ${blacklist} --numberOfProcessors ${task.cpus} --effectiveGenomeSize ${effective_genome_size}"
        if (params.reads_type == 'paired') {
            coverage_options += ' --extendReads'
        }
        if (params.assay == 'mnaseq') {
            coverage_options += ' --MNase'
        }
        if (params.assay == 'rnaseq') {
            // direction 'forward'/'reverse' keeps reads from that strand only;
            // anything else (e.g. 'unstranded') counts reads from both strands
            coverage_options = coverage_options.replace(' --extendReads', '')
            if (direction in ['forward', 'reverse']) {
                coverage_options += " --filterRNAstrand ${direction}"
            }
        }

        // RNA-seq makes one bigwig per strand from the same BAM, so the strand
        // goes into the name or the two files would overwrite each other.
        strand = (params.assay == 'rnaseq' && direction in ['forward', 'reverse']) ? "_${direction}" : ''
        if (scale_source == 'None') {    // unscaled track (a scaled one keeps its name even if the factor is 1)
            bw_name = "${bam.simpleName}${strand}_${params.normalize_by}.bw"
            coverage_options += " --scaleFactor 1 --normalizeUsing '${params.normalize_by}'"
        } else {
            bw_name = "${bam.simpleName}${strand}_${scale_source}-scaled.bw"
            coverage_options += " --normalizeUsing 'None' --scaleFactor ${scale_factor.toString().replace('\"', '')}"
        }

        """
        #!/bin/bash
        set -euo pipefail
        export MPLCONFIGDIR="\$PWD/matplotlib-config"
        mkdir -p "\$MPLCONFIGDIR"
        bamCoverage --bam ${bam} -o ${bw_name} ${coverage_options}
        """
}

process BIGWIG_BAMCOMPARE {
    label 'bamcompare'
    // retry if killed for memory, otherwise skip
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample_target.name} | vs ${bam_control.simpleName}"
    publishDir { "${readsDir(sample_target)}/bigwigs/" }, mode: params.publish_mode, overwrite: false, pattern: "*.bw"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control, stageAs: 'sample_control'), path(bam_control), path(bam_index_control)
        path(blacklist_file)
        val(effective_genome_size)
        path(genome_fasta)

    output:
        path "${bw_name}", emit: bigwig

    script:
        blacklist = (blacklist_file.name == 'NO_BLACKLIST_FILE') ? '' : "--blackListFileName ${blacklist_file}"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' ${blacklist} --numberOfProcessors ${task.cpus} --effectiveGenomeSize ${effective_genome_size}"
        if (params.reads_type == 'paired') {
            coverage_options += ' --extendReads'
        }
        if (params.assay == 'mnaseq') {
            coverage_options += ' --MNase'
        }
        bw_name = "${sample_target.simpleName}_without_${bam_control.simpleName}.bw"   // sample_control is staged as "sample_control", so name from the BAM

        """
        #!/bin/bash
        set -euo pipefail
        export MPLCONFIGDIR="\$PWD/matplotlib-config"
        mkdir -p "\$MPLCONFIGDIR"
        bamCompare -b1 ${bam_target} -b2 ${bam_control} -o ${bw_name} --operation subtract -of bigwig ${coverage_options} --effectiveGenomeSize ${effective_genome_size} --scaleFactorsMethod None
        bigWigToWig ${bw_name} ${bw_name.replaceAll('.bw', '.wig')}
        awk -v OFS='\t' -F'\t' '\$4 = (\$4 > 0 ? \$4 : 1e-30) 1' ${bw_name.replaceAll('.bw', '.wig')} > ${bw_name.replaceAll('.bw', '_nonegative.wig')}
        faSize ${genome_fasta} -detailed -tab > chrom.sizes
        wigToBigWig ${bw_name.replaceAll('.bw', '_nonegative.wig')} chrom.sizes ${bw_name}
        """
}

process BEDGRAPH_COVERAGE {
    label 'bedgraph'
    tag "${sample.name}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        path(blacklist_file)
        val(effective_genome_size)

    output:
        tuple path(sample), path("${bam.simpleName}.bedgraph"), emit: bedgraph

    script:
        blacklist = (blacklist_file.name == 'NO_BLACKLIST_FILE') ? '' : "--blackListFileName ${blacklist_file}"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' ${blacklist} --numberOfProcessors ${task.cpus} --effectiveGenomeSize ${effective_genome_size}"
        if (params.reads_type == 'paired') {
            coverage_options += ' --extendReads'
        }
        if (params.assay == 'mnaseq') {
            coverage_options += ' --MNase'
        }

        """
        #!/bin/bash
        set -euo pipefail
        export MPLCONFIGDIR="\$PWD/matplotlib-config"
        mkdir -p "\$MPLCONFIGDIR"
        bamCoverage --bam ${bam} --outFileFormat bedgraph -o ${bam.simpleName}.bedgraph ${coverage_options}
        """
}

process MACS_PEAKS {
    label 'macs_peaks'
    tag "${sample_target.name}"
    publishDir { "${readsDir(sample_target)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        val(effective_genome_size)

    output:
        tuple val(sample_target.name), path("${peaks_file}.narrowPeak"), emit: macs_narrowpeak
        tuple val(sample_target.name), path("${peaks_file}.broadPeak"), emit: macs_broadpeak
        tuple val(sample_target.name), path("${peaks_file}.gappedPeak"), emit: macs_gappedpeak

    script:
        bamformat  = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"
        parameters += macsShiftModelArgs(params.assay, params.reads_type)

        command    = "macs3 callpeak -t ${bam_target} ${parameters} --nolambda --outdir ./ -n ${bam_target.simpleName}"
        peaks_file = "${bam_target.simpleName}_peaks"
        """
        #!/bin/bash
        set -euo pipefail
        # narrow and broad calls run at the same time (each is single-threaded);
        # the broad call works in its own folder so the two don't share files
        rc=0
        ${command} > macs_narrow.log 2>&1 &
        narrow=\$!
        ${command} --broad --outdir broad > macs_broad.log 2>&1 &
        broad=\$!
        wait \$narrow || rc=\$?
        wait \$broad || rc=\$?
        if [ "\$rc" -ne 0 ]; then cat macs_narrow.log macs_broad.log >&2; exit "\$rc"; fi
        mv broad/${peaks_file}.broadPeak broad/${peaks_file}.gappedPeak .
        """
}

process MACS_PEAKS_WITH_CONTROL {
    label 'macs_peaks_control'
    tag "${sample_target.name} | vs ${bam_control.simpleName}"
    publishDir { "${readsDir(sample_target)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control, stageAs: 'sample_control'), path(bam_control), path(bam_index_control)
        val(effective_genome_size)

    output:
        tuple val(sample_target.name), path("${peaks_file}.narrowPeak"), emit: macs_narrowpeak
        tuple val(sample_target.name), path("${peaks_file}.broadPeak"), emit: macs_broadpeak
        tuple val(sample_target.name), path("${peaks_file}.gappedPeak"), emit: macs_gappedpeak

    script:
        bamformat  = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"
        parameters += macsShiftModelArgs(params.assay, params.reads_type)

        command    = "macs3 callpeak -t ${bam_target} -c ${bam_control} ${parameters} --outdir ./ -n ${bam_target.simpleName}_without_${bam_control.simpleName}-control"
        peaks_file = "${bam_target.simpleName}_without_${bam_control.simpleName}-control_peaks"
        """
        #!/bin/bash
        set -euo pipefail
        # narrow and broad calls run at the same time (each is single-threaded);
        # the broad call works in its own folder so the two don't share files
        rc=0
        ${command} > macs_narrow.log 2>&1 &
        narrow=\$!
        ${command} --broad --outdir broad > macs_broad.log 2>&1 &
        broad=\$!
        wait \$narrow || rc=\$?
        wait \$broad || rc=\$?
        if [ "\$rc" -ne 0 ]; then cat macs_narrow.log macs_broad.log >&2; exit "\$rc"; fi
        mv broad/${peaks_file}.broadPeak broad/${peaks_file}.gappedPeak .
        """
}

// Shared MACS3 --nomodel / --shift / --extsize logic used by both MACS_PEAKS processes
def macsShiftModelArgs(assay, reads_type) {
    if (assay == 'cutntag') {
        return (reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift 0 --extsize 200'
    } else if (assay == 'chipseq') {
        return ' --nomodel'
    } else {
        // atacseq and any other/default assay
        return (reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
    }
}

process SEACR_PEAKS {
    label 'seacr_peaks'
    // retry if killed for memory, otherwise skip
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample.name}"
    publishDir { "${readsDir(sample)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.stringent.bed"

    input:
        tuple path(sample), path(bedgraph)

    output:
        tuple val(sample.name), path("${peaks_file_prefix}.stringent.bed"), emit: seacr_peaks

    script:
        peaks_file_prefix = "${bedgraph.simpleName}_peaks"
        """
        #!/bin/bash
        set -euo pipefail
        SEACR ${bedgraph} ${params.seacr_threshold} non stringent ${peaks_file_prefix}
        """
}

process SEACR_PEAKS_WITH_CONTROL {
    label 'seacr_peaks_control'
    // retry if killed for memory, otherwise skip
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample_target.name} | vs ${bedgraph_control.simpleName}"
    publishDir { "${readsDir(sample_target)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.stringent.bed"

    input:
        tuple path(sample_target), path(bedgraph_target)
        tuple path(sample_control, stageAs: 'sample_control'), path(bedgraph_control)

    output:
        tuple val(sample_target.name), path("${peaks_file_prefix}.stringent.bed"), emit: seacr_peaks

    script:
        peaks_file_prefix = "${bedgraph_target.simpleName}_without_${bedgraph_control.simpleName}-control_peaks"
        """
        #!/bin/bash
        set -euo pipefail
        SEACR ${bedgraph_target} ${bedgraph_control} non stringent ${peaks_file_prefix}
        """
}

process GOPEAKS_PEAKS {
    label 'gopeaks_peaks'
    // retry if killed for memory, otherwise skip (like the other secondary peak callers)
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample.name}"
    publishDir { "${readsDir(sample)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.bed"

    input:
        tuple path(sample), path(bam), path(bam_index)

    output:
        tuple val(sample.name), path("${bam.simpleName}_gopeaks_peaks.bed"), emit: gopeaks_peaks

    script:
        """
        #!/bin/bash
        set -euo pipefail
        gopeaks -b ${bam} -o ${bam.simpleName}_gopeaks -p ${params.gopeaks_pvalue}
        """
}

process GOPEAKS_PEAKS_WITH_CONTROL {
    label 'gopeaks_peaks_control'
    // retry if killed for memory, otherwise skip (like the other secondary peak callers)
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample_target.name} | vs ${bam_control.simpleName}"
    publishDir { "${readsDir(sample_target)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.bed"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control, stageAs: 'sample_control'), path(bam_control), path(bam_index_control)

    output:
        tuple val(sample_target.name), path("${peaks_file}_gopeaks_peaks.bed"), emit: gopeaks_peaks

    script:
        peaks_file = "${bam_target.simpleName}_without_${bam_control.simpleName}-control"
        """
        #!/bin/bash
        set -euo pipefail
        gopeaks -b ${bam_target} -c ${bam_control} -o ${bam_target.simpleName}_without_${bam_control.simpleName}-control_gopeaks -p ${params.gopeaks_pvalue}
        """
}

// Genrich reads alignments grouped by read name, so each final BAM is sorted
// by name once here (multi-threaded) and shared by both Genrich processes.
process NAME_SORT_BAM {
    label 'name_sort_bam'
    tag "${sample.name} | ${bam.name}"

    input:
        tuple path(sample), path(bam), path(bam_index)

    output:
        tuple path(sample), path("${bam.simpleName}.namesorted.bam"), emit: bam_namesorted

    script:
        """
        #!/bin/bash
        set -euo pipefail
        samtools sort -n -@ ${task.cpus} -o ${bam.simpleName}.namesorted.bam ${bam}
        """
}

// Genrich options shared by both Genrich processes, chosen from the assay and
// read type:
//  - ATAC-seq: -j (ATAC mode: intervals centred on Tn5 cut sites, with the
//    standard +5/-5 bp Tn5 shift); single-end reads also need -y so unpaired
//    alignments (their 5' cut site) are used.
//  - Other assays, paired-end: whole fragments from proper pairs (default).
//  - Other assays, single-end: -w extends reads to the fragment size, like
//    MACS --extsize.
// BAMs are already MAPQ-filtered and de-duplicated, so -m and -r are not used.
// Text value of an optional param, '' when unset. An empty value on the
// command line (e.g. --genrich_exclude_chr '') arrives as boolean true.
def optionalParam(value) {
    def text = (value == null || value instanceof Boolean) ? '' : value.toString().trim()
    return (text in ['true', 'false']) ? '' : text
}

def genrichArgs(blacklist_file) {
    def args = []
    if (params.assay == 'atacseq') {
        args << '-j'
        if (params.reads_type != 'paired') { args << '-y' }
    } else if (params.reads_type != 'paired') {
        args << "-w ${params.genrich_extsize}"
    }
    def exclude_chr = optionalParam(params.genrich_exclude_chr)
    def qvalue      = optionalParam(params.genrich_qvalue)
    if (exclude_chr) { args << "-e ${exclude_chr}" }
    if (blacklist_file.name != 'NO_BLACKLIST_FILE') { args << "-E ${blacklist_file}" }
    args << (qvalue ? "-q ${qvalue}" : "-p ${params.genrich_pvalue}")
    args << "-a ${params.genrich_min_auc}"
    args << "-g ${params.genrich_max_gap}"
    return args.join(' ')
}

// Run Genrich with its -v progress/count report going to a log for QC. If it
// fails, the log is printed to stderr and its real exit code kept, so an
// out-of-memory kill is retried with more memory.
def genrichCommand(genrich_args, log_file) {
    return [
        "set +e",
        "Genrich ${genrich_args} -v 2> ${log_file}",
        'rc=$?',
        "set -e",
        'if [ "$rc" -ne 0 ]; then',
        "    cat ${log_file} >&2",
        '    exit "$rc"',
        "fi"
    ].join('\n        ')
}

process GENRICH_PEAKS {
    label 'genrich_peaks'
    // retry if killed for memory, otherwise skip (like the other secondary peak callers)
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample.name}"
    publishDir { "${readsDir(sample)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.narrowPeak"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "genrich_*.log"

    input:
        tuple path(sample), path(bam)
        path(blacklist_file)

    output:
        tuple val(sample.name), path("${peaks_name}_genrich_peaks.narrowPeak"), emit: genrich_peaks
        path "genrich_${peaks_name}.log"

    script:
        peaks_name = "${bam.simpleName}"
        genrich_cmd = genrichCommand("-t ${bam} -o ${peaks_name}_genrich_peaks.narrowPeak ${genrichArgs(blacklist_file)}", "genrich_${peaks_name}.log")
        """
        #!/bin/bash
        set -euo pipefail
        ${genrich_cmd}
        """
}

process GENRICH_PEAKS_WITH_CONTROL {
    label 'genrich_peaks_control'
    errorStrategy { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'ignore' }
    tag "${sample_target.name} | vs ${bam_control.simpleName}"
    publishDir { "${readsDir(sample_target)}/peaks/" }, mode: params.publish_mode, overwrite: false, pattern: "*.narrowPeak"
    publishDir { "${readsDir(sample_target)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "genrich_*.log"

    input:
        tuple path(sample_target), path(bam_target)
        tuple path(sample_control, stageAs: 'sample_control'), path(bam_control)
        path(blacklist_file)

    output:
        tuple val(sample_target.name), path("${peaks_name}_genrich_peaks.narrowPeak"), emit: genrich_peaks
        path "genrich_${peaks_name}.log"

    script:
        peaks_name = "${bam_target.simpleName}_without_${bam_control.simpleName}-control"
        genrich_cmd = genrichCommand("-t ${bam_target} -c ${bam_control} -o ${peaks_name}_genrich_peaks.narrowPeak ${genrichArgs(blacklist_file)}", "genrich_${peaks_name}.log")
        """
        #!/bin/bash
        set -euo pipefail
        ${genrich_cmd}
        """
}

// FRiP (fraction of reads in peaks) for every peak file a sample produced,
// scored against that sample's final BAM, written to one table in QC/.
// Reads = mapped primary alignments (-F 0x904: no unmapped, secondary or
// supplementary); paired-end mates are counted individually, as in ENCODE's
// FRiP. A read overlapping several peaks is counted once.
process FRIP {
    label 'frip'
    tag "${sample.name}"
    publishDir { "${readsDir(sample)}/QC/" }, mode: params.publish_mode, overwrite: false, pattern: "*_FRiP.txt"

    input:
        tuple path(sample), path(bam), path(bam_index), path(peak_files, stageAs: 'peaks/*')

    output:
        path "${sample.simpleName}_FRiP.txt"

    script:
        """
        #!/bin/bash
        set -euo pipefail
        out=${sample.simpleName}_FRiP.txt
        total=\$(samtools view -c -@ ${task.cpus} -F 0x904 ${bam})
        {
            echo "# FRiP (fraction of reads in peaks) for ${sample.name}"
            echo "# Reads: ${bam}, mapped primary alignments (samtools -F 0x904); paired-end mates counted individually."
            echo "# A read overlapping several peaks counts once. total_reads includes reads on chromosomes a caller excluded (e.g. chrM)."
            printf 'peak_file\\tpeak_caller\\tpeaks\\treads_in_peaks\\ttotal_reads\\tFRiP\\n'
            for f in \$(ls peaks/ | sort); do
                case "\$f" in
                    *_genrich_peaks.narrowPeak) caller=Genrich ;;
                    *_gopeaks_peaks.bed)        caller=GoPeaks ;;
                    *.stringent.bed)            caller=SEACR ;;
                    *.narrowPeak)               caller=MACS3_narrow ;;
                    *.broadPeak)                caller=MACS3_broad ;;
                    *.gappedPeak)               caller=MACS3_gapped ;;
                    *)                          caller=other ;;
                esac
                peaks=\$(awk '!/^(#|track|browser)/ && NF >= 3 {n++} END {print n+0}' "peaks/\$f")
                in_peaks=\$(samtools view -c -@ ${task.cpus} -F 0x904 -L "peaks/\$f" ${bam})
                frip=\$(awk -v a="\$in_peaks" -v b="\$total" 'BEGIN { if (b > 0) printf "%.4f", a / b; else print "NA" }')
                printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "\$f" "\$caller" "\$peaks" "\$in_peaks" "\$total" "\$frip"
            done
        } > "\$out"
        """
}
