/*
    -------------------------------------------------
        VERIFY AND VALIDATE READS FILES FOR INPUT
    -------------------------------------------------
*/
process GET_READS {
    input:
        val reads_dir
    output:
        path reads, emit: raw_reads
    script:
        println "${reads_dir.getProperties()}"
        if ("${reads_dir}" != ''){
            assert file("${reads_dir}").exists() : "Cannot find target reads directory ${reads_dir}"
            if (params.reads_type == 'paired'){
                assert files("${reads_dir}/*_R{1,2}.fastq.gz").size() > 0 : "Expected paired reads...check files in ${reads_dir}"
                assert files("${reads_dir}/*_R{1,2}.fastq.gz").size() %2 == 0 : "Expected paired reads...check files in ${reads_dir}"
                reads = Channel.fromPath("${reads_dir}/*_R{1,2}.fastq.gz", type: 'file')
            } else {
                assert files("${reads_dir}/*.fastq.gz").size() > 0 : "Expected reads...check files in ${reads_dir}"
                reads = Channel.fromPath("${reads_dir}/*.fastq.gz", type: 'file')
            }
        } else {
            reads = Channel.empty()
        }
}

/*
    -------------------------------------
        DEFINE PROCESSES FOR PIPELINE
    -------------------------------------
*/
// Validate integrity of reads files
process MD5SUMCHECK {
    cpus 1
    label 'md5sum_check'
    tag "md5sum_check_${read}"
    publishDir "${sample}", mode: 'copy', overwrite: false, pattern: "*.{md5}"

    input:
        tuple path(sample), path(read), path(md5)

    output:
        file md5 

    script:
        if (!md5.exists()) {
            """
            #!/bin/sh
            set -e
            echo "Generating md5 checksum hash for ${read}..."
            md5sum ${read} > ${md5};
            md5sum -c ${md5}
            """
        } else {
            """
            #!/bin/sh
            set -e
            echo "Checking md5 checksum hash for ${read}..."
            md5sum -c ${md5}
            """
        }
}

process FASTQC {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.2 as int) 
    label 'fastqc_check'
    tag "fastqc_check_${sample}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "*fastqc.{zip,html}"

    input:
        tuple path(sample), path(reads)

    output:
        path "*fastqc.{zip,html}"

    script:
        """
        #!/bin/bash
        fastqc -t $task.cpus ${reads}
        """
}

process CUTADAPT_PAIRED {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.25 as int) 
    label 'cutadapt_trim_paired'
    tag "cutadapt_trim_paired_${sample.baseName}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "cutadapt_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        path "cutadapt_${sample.baseName}.log"
        //tuple path(sample), path("${reads.first.name.replace('.fastq.gz', '_trimmed.fastq')}"), path("${reads.last.name.replace('.fastq.gz', '_trimmed.fastq')}"), emit: trimmed_reads
        tuple path(sample), path("*_R{1,2}_trimmed.fastq"), emit: trimmed_reads

    script:
        adapters = params.adapters ? "-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT" : ""
        trim = params.qctrim ? "-m 20 -q 20,20 -u 11 -U 11" : ""

        if (params.assay == 'rnaseq' & params.polyAtrim) {
            """
            #!/bin/bash
            cutadapt --cores=${task.cpus} $trim $adapters -o ${reads[0].name.replace('.fastq.gz', '_adapt-trimmed.fastq')} -p ${reads[1].name.replace('.fastq.gz', '_adapt-trimmed.fastq')} ${reads[0]} ${reads[1]} &>> cutadapt_${sample.baseName}.log
            cutadapt --cores=${task.cpus} --poly-a -o ${reads[0].name.replace('.fastq.gz', '_trimmed.fastq')} -p ${reads[1].name.replace('.fastq.gz', '_trimmed.fastq')} ${reads[0].name.replace('.fastq.gz', '_adapt-trimmed.fastq')} ${reads.name[1].replace('.fastq.gz', '_adapt-trimmed.fastq')} &>> cutadapt_${sample.baseName}.log
            """
        } else {
            """
            #!/bin/bash
            cutadapt --cores=${task.cpus} $trim $adapters -o ${reads[0].name.replace('.fastq.gz', '_trimmed.fastq')} -p ${reads[1].name.replace('.fastq.gz', '_trimmed.fastq')} ${reads[0]} ${reads[1]} &>> cutadapt_${sample.baseName}.log
            """
        }

}
// TODO: MODIFY FILE NAMES TO THE INPUT FILE NAME INSTEAD OF THE SAMPLE DIRECTORY NAME AS IN THIS PROCESS...
process CUTADAPT_SINGLE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.2 as int) 
    label 'cutadapt_trim_single'
    tag "cutadapt_trim_single_${reads.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "cutadapt_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        //path "${reads.name.replaceAll('.fastq.gz', '_trimmed.fastq')}"
        path "cutadapt_${reads.name.split("\\.")[0]}.log"
        tuple path(sample), path("${reads.name.replace('.fastq.gz', '_trimmed.fastq')}"), emit: trimmed_reads

    script:
        adapters = params.adapters ? '-a CTGTCTCTTATACACATCT' : ''
        trim = params.qctrim ? '-m 20 -q 20 -u 11' : ''

        if (params.assay == 'rnaseq' & params.polyAtrim) {
            """
            #!/bin/bash
            cutadapt --cores=${task.cpus} $trim $adapters -o ${reads.name.replace('.fastq.gz', '_adapt-trimmed.fastq')} $reads &>> cutadapt_${reads.name.split("\\.")[0]}.log
            cutadapt --cores=${task.cpus} --poly-a -o ${reads.name.replace('.fastq.gz', '_trimmed.fastq')} ${reads.name.replace('.fastq.gz', '_adapt-trimmed.fastq')} &>> cutadapt_${reads.name.split("\\.")[0]}.log
            """
        } else {
            """
            #!/bin/bash
            cutadapt --cores=${task.cpus} $trim $adapters -o ${reads.name.replace('.fastq.gz', '_trimmed.fastq')} $reads &>> cutadapt_${reads.name.split("\\.")[0]}.log
            """
        }
}


process BOWTIE2_PAIRED {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bowtie2_align_paired'
    tag "bowtie2_align_paired_${sample.baseName}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "bowtie2_*.log"

    input:
        tuple path(sample), path(reads)
    
    output:
        path "bowtie2_${sample.baseName}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads
    
    script:
        reads_name = ["${reads[0].name}", "${reads[1].name}"]
            .collect { it.replaceAll('_R[1|2]_trimmed.fastq', '').replaceAll('_R[1|2].fastq', '')}
            .unique().join()
        """
        #!/bin/bash
        bowtie2 -p $task.cpus --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x $params.genome_index -1 $reads[0] -2 $reads[1] \
        --rg-id $params.rg_id --rg SM:$params.rg_sm --rg LB:$params.rg_lb --rg PU:$params.rg_pu --rg PL:$params.rg_pl \
        2> bowtie2_${sample.baseName}.log | samtools view -bS - > ${reads_name}.bam 
        """
}

process BOWTIE2_SINGLE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bowtie2_align_single'
    tag "bowtie2_align_single_${reads.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "bowtie2_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        path "bowtie2_${reads.name.split("\\.")[0]}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads

    script:
        reads_name = reads.name.replace('_trimmed.fastq', '').replace('.fastq', '')
        """
        #!/bin/bash
        bowtie2 -p $task.cpus --local --very-sensitive-local --no-unal --phred33 \
        -x $params.genome_index -U $reads \
        --rg-id $params.rg_id --rg SM:$params.rg_sm --rg LB:$params.rg_lb --rg PU:$params.rg_pu --rg PL:$params.rg_pl \
        2> bowtie2_${reads.name.split("\\.")[0]}.log | samtools view -bS - > ${reads_name}.bam
        """
}


process HISAT2_PAIRED {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'hisat2_align_paired'
    tag "hisat2_align_paired_${sample.baseName}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "hisat2_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        path "hisat2_${sample.baseName}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads

    script:
        reads_name = ["${reads[0].name}", "${reads[1].name}"]
            .collect { it.replaceAll('_R[1|2]_trimmed.fastq', '').replaceAll('_R[1|2].fastq', '')}
            .unique().join()
        """
        #!/bin/bash
        hisat2 -p $task.cpus --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -x $params.genome_index -1 $reads[0] -2 $reads[1] \
        --rg-id $params.rg_id --rg SM:$params.rg_sm --rg LB:$params.rg_lb --rg PU:$params.rg_pu --rg PL:$params.rg_pl \
        2> hisat2_${sample.baseName}.log | samtools view -bS - > ${reads_name}.bam
        """
}

process HISAT2_SINGLE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'hisat2_align_single'
    tag "hisat2_align_single_${reads.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "hisat2_*.log"

    input:
        tuple path(sample), path(reads)

    output:
        path "hisat2_${reads.name.split("\\.")[0]}.log"
        tuple path(sample), path("${reads_name}.bam"), emit: mapped_reads

    script:
        reads_name = reads.name.replace('_trimmed.fastq', '').replace('.fastq', '')
        """
        #!/bin/bash
        hisat2 -p $task.cpus --no-unal --phred33 \
        -x $params.genome_index -U $reads \
        --rg-id $params.rg_id --rg SM:$params.rg_sm --rg LB:$params.rg_lb --rg PU:$params.rg_pu --rg PL:$params.rg_pl \
        2> hisat2_${reads.name.split("\\.")[0]}.log | samtools view -bS - > ${reads_name}.bam
        """
}


process SORT_BAM {
    cpus 1
    label 'sort_bam'
    tag "sort_bam_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "*.log"

    input:
        tuple path(sample), path(bam)

    output:
        path "picard-sort_${bam.name.split("\\.")[0]}.log"
        path "picard-dupStats_${bam.name.split("\\.")[0]}.log"
        path "samtools-flagstat-prefiltered_${bam.name.split("\\.")[0]}.log"
        tuple path(sample), path("${bam.name.replace('.bam', '.sorted.bam')}"), emit: bam_sorted

    script:
        """
        #!/bin/bash
        picard SortSam -I $bam -O ${bam.name.replace('.bam', '.sorted.bam')} -SORT_ORDER $params.sort_bam_order &>> picard-sort_${bam.name.split("\\.")[0]}.log
        samtools flagstat ${bam.name.replace('.bam', '.sorted.bam')} &>> samtools-flagstat-prefiltered_${bam.name.split("\\.")[0]}.log
        picard MarkDuplicates -I ${bam.name.replace('.bam', '.sorted.bam')} -O ${bam.name.replace('.bam', '.dupMarked.bam')} -METRICS_FILE picard-dupStats_${bam.name.split("\\.")[0]}.log
        """
}

process QFILTER_BAM {
    cpus 1
    label 'qfilter_bam'
    tag "qfilter_bam_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "*.log"

    input:
        tuple path(sample), path(bam)

    output:
        path "samtools-sort_${bam.name.split("\\.")[0]}.log"
        path "samtools-flagstat-postqfiltered_${bam.name.split("\\.")[0]}.log"
        tuple path(sample), path("${bam.name.replace('.sorted.bam', '')}.MAPQ${params.mapq}.bam"), emit: bam_qfiltered

    script:
        sam_flags = params.reads_type == 'paired' ? '-f 2' : '-F 4'
        """
        #!/bin/bash
        samtools view -bu $sam_flags $bam | samtools view -b -q $params.mapq - | samtools sort - -o ${bam.name.replace('.sorted.bam', '')}.MAPQ${params.mapq}.bam &>> samtools-sort_${bam.name.split("\\.")[0]}.log
        samtools flagstat ${bam.name.replace('.sorted.bam', '')}.MAPQ${params.mapq}.bam &>> samtools-flagstat-postqfiltered_${bam.name.split("\\.")[0]}.log
        """
}

process DEDUPLICATE_BAM {
    cpus 1
    label 'deduplicate_bam'
    tag "deduplicate_bam_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/QC/", mode: 'copy', overwrite: false, pattern: "*.log"

    input:
        tuple path(sample), path(bam)

    output:
        path "picard-deduplicate_${bam.name.split("\\.")[0]}.log"
        path "samtools-flagstat_${bam.name.split("\\.")[0]}.log"
        tuple path(sample), path("${bam.name.replace('.bam', '.NoDups.bam')}"), emit: bam_deduplicated

    script:
        """
        #!/bin/bash
        picard MarkDuplicates -I ${bam} -O ${bam.name.replace('.bam', '.NoDups.bam')} -REMOVE_DUPLICATES true -METRICS_FILE picard-deduplicate_${bam.name.split("\\.")[0]}.log
        samtools flagstat ${bam.name.replace('.bam', '.NoDups.bam')} &>> samtools-flagstat_${bam.name.split("\\.")[0]}.log
        """
}

process INDEX_BAM {
    cpus 1
    label 'index_bam'
    tag "index_bam_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/bams/", mode: 'copy', overwrite: false, pattern: "*.{bam,bai}"

    input:
        tuple path(sample), path(bam)

    output:
        tuple path(sample), path("${bam}"), path("${bam}.bai"), emit: bam_indexed

    script:
        """
        #!/bin/bash
        samtools index $bam
        """
}

process MERGE_BAMS {
    cpus 1
    label 'merge_bams'
    tag "merge_bams_${sample.baseName}"
    publishDir "${sample}/bams/", mode: 'copy', overwrite: false, pattern: "*{.bam,output.txt}"

    input:
        tuple path(sample), path(bams), path(bam_index)

    output:
        path "picard-merge_${sample.baseName}.log"
        tuple path(sample), path("${sample.baseName}.bam"), path("${sample.baseName}.bam.bai"), emit: bam_merged

    script:
        """
        #!/bin/bash
        picard MergeSamFiles $bams O=${sample.baseName}.bam &>> picard-merge_${sample.baseName}.log
        samtools index $bams
        """
}


process MULTIQC {
    cpus 1
    label 'multiqc_report'
    tag "multiqc_report_${sample.baseName}"
    publishDir "${sample}/", mode: 'copy', overwrite: false, pattern: "*.{html}"

    input:
        tuple path(sample), path(bam), path(bai)
        path multiqc_config

    output:
        path "${sample.baseName}_multiqc_report.html"

    script:
        config_file = ("$multiqc_config" == "NO_MULTIQC_CONFIG_FILE") ? "" : "-c $multiqc_config"
        """
        #!/bin/bash
        multiqc ${sample}/ $config_file -v --force --filename ${sample.baseName}_multiqc_report.html
        """
}

process BIGWIG_COVERAGE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bigwig'
    tag "bigwig_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/bigwigs/", mode: 'copy', overwrite: false, pattern: "*.{bw}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        path blacklist_file
        val effective_genome_size

    output:
        path "${bam.name.split('\\.')[0]}_${params.normalize_by}.bw", emit: bigwig

    script:
        blacklist = ("$blacklist_file" == "NO_BLACKLIST_FILE") ? "" : "--blackListFileName $blacklist_file"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' $blacklist --normalizeUsing '${params.normalize_by}' --numberOfProcessors $task.cpus"
        if (params.reads_type == 'paired'){
            coverage_options += " --extendReads"
        }
        if (params.assay == 'mnaseq'){
            coverage_options += ' --MNase'
        }
        if (params.assay == 'rnaseq'){
            coverage_options = coverage_options.replace(" --extendReads", "")
        }
        
        """
        #!/bin/bash
        bamCoverage --bam $bam -o ${bam.name.split('\\.')[0]}_${params.normalize_by}.bw $coverage_options --effectiveGenomeSize $effective_genome_size
        """
}

process BIGWIG_COVERAGE_STRANDED {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bigwig_stranded'
    tag "bigwig_stranded_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/bigwigs/", mode: 'copy', overwrite: false, pattern: "*.{bw}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        path blacklist_file
        val effective_genome_size
        val direction

    output:
        path "${bam.name.split('\\.')[0]}_${params.normalize_by}.${direction}.bw", emit: bigwig

    script:
        blacklist = ("$blacklist_file" == "NO_BLACKLIST_FILE") ? "" : "--blackListFileName $blacklist_file"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' $blacklist --normalizeUsing '${params.normalize_by}' --numberOfProcessors $task.cpus"
        """
        #!/bin/bash
        bamCoverage --bam $bam -o ${bam.name.split('\\.')[0]}_${params.normalize_by}.${direction}.bw $coverage_options --effectiveGenomeSize $effective_genome_size --filterRNAstrand $direction
        """
}

process BIGWIG_BAMCOMPARE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bamcompare'
    tag "bamcompare_${sample_target.baseName}"
    publishDir "${sample_target}/bigwigs/", mode: 'copy', overwrite: false, pattern: "*.{bw}"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control), path(bam_control), path(bam_index_control)
        path blacklist_file
        val effective_genome_size

    output:
        path "${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.bw", emit: bigwig

    script:
        blacklist = ("$blacklist_file" == "NO_BLACKLIST_FILE") ? "" : "--blackListFileName $blacklist_file"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' $blacklist --normalizeUsing '${params.normalize_by}' --numberOfProcessors $task.cpus"
        if (params.reads_type == 'paired'){
            coverage_options += " --extendReads"
        }
        if (params.assay == 'mnaseq'){
            coverage_options += ' --MNase'
        }
        if (params.assay == 'rnaseq'){
            coverage_options = coverage_options.replace(" --extendReads", "")
        }
        """
        #!/bin/bash
        bamCompare -b1 $bam_target -b2 $bam_control -o ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.bw --operation subtract -of bigwig $coverage_options --effectiveGenomeSize $effective_genome_size --scaleFactorsMethod None
        bigWigToWig ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.bw ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.wig
        awk -v OFS='\t' -F'\t' '\$4 = (\$4 > 0 ? \$4 : 1e-30) 1' ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.wig > ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}_nonegative.wig
        faSize ${params.genome_index}.fa -detailed -tab > chrom.sizes
        wigToBigWig ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}_nonegative.wig chrom.sizes ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}_${params.normalize_by}.bw
        """
}

process BEDGRAPH_COVERAGE {
    cpus Math.max(1, Runtime.runtime.availableProcessors() * 0.5 as int) 
    label 'bedgraph'
    tag "bedgraph_${bam.name.split("\\.")[0]}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        path blacklist_file
        val effective_genome_size

    output:
        tuple path(sample), path("${bam.name.split('\\.')[0]}_${params.normalize_by}.bedgraph"), emit: bedgraph

    script:
        blacklist = ("$blacklist_file" == "NO_BLACKLIST_FILE") ? "" : "--blackListFileName $blacklist_file"
        coverage_options = "--binSize 10 --ignoreForNormalization 'chrM' $blacklist --normalizeUsing '${params.normalize_by}' --numberOfProcessors $task.cpus"
        if (params.reads_type == 'paired'){
            coverage_options += " --extendReads"
        }
        if (params.assay == 'mnaseq'){
            coverage_options += ' --MNase'
        }
        if (params.assay == 'rnaseq'){
            coverage_options = coverage_options.replace(" --extendReads", "")
        }
        
        """
        #!/bin/bash
        bamCoverage --bam $bam --outFileFormat bedgraph -o ${bam.name.split('\\.')[0]}_${params.normalize_by}.bedgraph $coverage_options --effectiveGenomeSize $effective_genome_size
        """
        
}


// CONSIDER SPIKEIN/NOSPIKEIN...?
// review https://macs3-project.github.io/MACS/docs/Advanced_Step-by-step_Peak_Calling.html

process MACS_NARROWPEAKS_NO_CONTROL {
    cpus 1
    label 'macs_narrowpeak_calling'
    tag "macs_narrowpeak_calling_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{narrowPeak}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        val (effective_genome_size)

    output:
        path "${peaks_file}.narrowPeak", emit: macs_narrowpeak

    script:
        bamformat = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"

        if (params.assay == 'cutntag'){
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift 0 --extsize 200'
        } else if (params.assay == 'chipseq') {
            parameters += " --nomodel"
        } else if (params.assay == 'atacseq') {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        } else {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        }

        command = "macs3 callpeak -t ${bam} $parameters --nolambda --outdir ./ -n ${bam.name.split('\\.')[0]}"
        peaks_file = "${bam.name.split('\\.')[0]}_peaks"
        """
        #!/bin/bash
        $command
        """
}

process MACS_NARROWPEAKS_CONTROL {
    cpus 1
    label 'macs_narrowpeak_calling'
    tag "macs_narrowpeak_calling_${sample_target.baseName}"
    publishDir "${sample_target}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{narrowPeak}"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control), path(bam_control), path(bam_index_control)
        val (effective_genome_size)

    output:
        path "${peaks_file}.narrowPeak", emit: macs_narrowpeak

    script:
        bamformat = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"

        if (params.assay == 'cutntag'){
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift 0 --extsize 200'
        } else if (params.assay == 'chipseq') {
            parameters += " --nomodel"
        } else if (params.assay == 'atacseq') {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        } else {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        }

        command = "macs3 callpeak -t ${bam_target} -c ${bam_control} $parameters --outdir ./ -n ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control"
        peaks_file = "${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control_peaks"
        
        """
        #!/bin/bash
        $command
        """
}

process MACS_BROADPEAKS_NO_CONTROL {
    cpus 1
    label 'macs_broadpeak_calling'
    tag "macs_broadpeak_calling_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{broadPeak}"

    input:
        tuple path(sample), path(bam), path(bam_index)
        val (effective_genome_size)

    output:
        path "${peaks_file}.broadPeak", emit: macs_broadpeak

    script:
        bamformat = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"

        if (params.assay == 'cutntag'){
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift 0 --extsize 200'
        } else if (params.assay == 'chipseq') {
            parameters += " --nomodel"
        } else if (params.assay == 'atacseq') {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        } else {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        }

        command = "macs3 callpeak -t ${bam} $parameters --nolambda --outdir ./ -n ${bam.name.split('\\.')[0]}"
        peaks_file = "${bam.name.split('\\.')[0]}_peaks"
        """
        #!/bin/bash
        $command --broad
        """
}

process MACS_BROADPEAKS_CONTROL {
    cpus 1
    label 'macs_broadpeak_calling'
    tag "macs_broadpeak_calling_${sample_target.baseName}"
    publishDir "${sample_target}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{broadPeak}"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control), path(bam_control), path(bam_index_control)
        val (effective_genome_size)

    output:
        path "${peaks_file}.broadPeak", emit: macs_broadpeak

    script:
        bamformat = (params.reads_type == 'paired') ? 'BAMPE' : 'BAM'
        parameters = "--keep-dup all --format ${bamformat} --gsize ${effective_genome_size} --qvalue ${params.macs_qvalue}"

        if (params.assay == 'cutntag'){
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift 0 --extsize 200'
        } else if (params.assay == 'chipseq') {
            parameters += " --nomodel"
        } else if (params.assay == 'atacseq') {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        } else {
            parameters += (params.reads_type == 'paired') ? ' --nomodel' : ' --nomodel --shift -100 --extsize 200'
        }

        command = "macs3 callpeak -t ${bam_target} -c ${bam_control} $parameters --outdir ./ -n ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control"
        peaks_file = "${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control_peaks"
        
        """
        #!/bin/bash
        $command --broad
        """
}

process SEACR_PEAKS_CONTROL {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3
    label 'seacr_peak_calling'
    tag "seacr_peak_calling_${sample_target.baseName}"
    publishDir "${sample_target}/peaks/", mode: 'copy', overwrite: false, pattern: "*.stringent.{bed}"

    input:
        tuple path(sample_target), path(bedgraph_target)
        tuple path(sample_control), path(bedgraph_control)

    output:
        path "${peaks_file_prefix}.stringent.bed", emit: seacr_peaks

    script:
        peaks_file_prefix = "${bedgraph_target.name.split('\\.')[0]}_without_${bedgraph_control.name.split('\\.')[0]}-control_peaks"
        """
        #!/bin/bash
        SEACR_1.3.sh ${bedgraph_target} ${bedgraph_control} non stringent ${peaks_file_prefix}
        """
}

process SEACR_PEAKS_NO_CONTROL {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3
    label 'seacr_peak_calling'
    tag "seacr_peak_calling_${bedgraph.name.split("\\.")[0]}"
    publishDir "${sample}/peaks/", mode: 'copy', overwrite: false, pattern: "*.stringent.{bed}"

    input:
        tuple path(sample), path(bedgraph)

    output:
        path "${peaks_file_prefix}.stringent.bed", emit: seacr_peaks

    script:
        peaks_file_prefix = "${bedgraph.name.split('\\.')[0]}_peaks"
        """
        #!/bin/bash
        SEACR_1.3.sh ${bedgraph} ${params.seacr_threshold} non stringent ${peaks_file_prefix}
        """
}

process GOPEAKS_PEAKS_CONTROL {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3
    memory { 10.GB * task.attempt }
    label 'gopeaks_peak_calling'
    tag "gopeaks_peak_calling_${sample_target.baseName}"
    publishDir "${sample_target}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{bed}"

    input:
        tuple path(sample_target), path(bam_target), path(bam_index_target)
        tuple path(sample_control), path(bam_control), path(bam_index_control)

    output:
        path "${peaks_file}_gopeaks_peaks.bed", emit: gopeaks_peaks

    script:
        peaks_file = "${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control"
        """
        #!/bin/bash
        gopeaks -b ${bam_target} -c ${bam_control} -o ${bam_target.name.split('\\.')[0]}_without_${bam_control.name.split('\\.')[0]}-control_gopeaks -p ${params.gopeaks_pvalue}
        """
}

process GOPEAKS_PEAKS_NO_CONTROL {
    cpus 1
    errorStrategy 'retry'
    maxRetries 3
    memory { 10.GB * task.attempt }
    label 'gopeaks_peak_calling'
    tag "gopeaks_peak_calling_${bam.name.split("\\.")[0]}"
    publishDir "${sample}/peaks/", mode: 'copy', overwrite: false, pattern: "*.{bed}"

    input:
        tuple path(sample), path(bam), path(bam_index)

    output:
        path "${peaks_file}_gopeaks_peaks.bed", emit: gopeaks_peaks

    script:
        peaks_file = "${bam.name.split('\\.')[0]}"
        """
        #!/bin/bash
        gopeaks -b ${bam} -o ${bam.name.split('\\.')[0]}_gopeaks -p ${params.gopeaks_pvalue}
        """
}