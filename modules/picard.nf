#!/usr/bin/env nextflow

process COLLECT_QUALITY_YIELD_METRICS {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory '3.5 GB'

    input:
    tuple val(output_prefix), path(bam)

    output:
    path "${output_prefix}.quality_yield_metrics.txt", emit: quality_yield_metrics

    script:
    """
    java -Xms2000m -Xmx3000m -jar /usr/picard/picard.jar \
        CollectQualityYieldMetrics \
        --INPUT ${bam} \
        --USE_ORIGINAL_QUALITIES true \
        --OUTPUT ${output_prefix}.quality_yield_metrics.txt
    """
}

process COLLECT_READ_GROUP_BAM_METRICS {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 7.GB

    input:
    tuple val(output_prefix), path(bam)
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    val collect_gc_bias_metrics

    output:
    path "${output_prefix}.alignment_summary_metrics", emit: alignment_summary_metrics
    path "${output_prefix}.gc_bias.detail_metrics", emit: gc_bias_detail_metrics, optional: true
    path "${output_prefix}.gc_bias.pdf", emit: gc_bias_pdf, optional: true
    path "${output_prefix}.gc_bias.summary_metrics", emit: gc_bias_summary_metrics, optional: true

    script:
    def collect_gc_bias_metrics_arg = collect_gc_bias_metrics ? "--PROGRAM CollectGcBiasMetrics" : ""
    """
    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
        CollectMultipleMetrics \
        --INPUT ${bam} \
        --REFERENCE_SEQUENCE ${ref_fasta} \
        --OUTPUT ${output_prefix} \
        --ASSUME_SORTED true \
        --PROGRAM null \
        --PROGRAM CollectAlignmentSummaryMetrics \
        ${collect_gc_bias_metrics_arg} \
        --METRIC_ACCUMULATION_LEVEL null \
        --METRIC_ACCUMULATION_LEVEL READ_GROUP
    """
}

process COLLECT_UNSORTED_READ_GROUP_BAM_METRICS {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 7.GB

    input:
    tuple val(output_prefix), path(bam)

    output:
    path "${output_prefix}.base_distribution_by_cycle.pdf", emit: base_distribution_by_cycle_pdf
    path "${output_prefix}.base_distribution_by_cycle_metrics", emit: base_distribution_by_cycle_metrics
    path "${output_prefix}.insert_size_histogram.pdf", emit: insert_size_histogram_pdf
    path "${output_prefix}.insert_size_metrics", emit: insert_size_metrics
    path "${output_prefix}.quality_by_cycle.pdf", emit: quality_by_cycle_pdf
    path "${output_prefix}.quality_by_cycle_metrics", emit: quality_by_cycle_metrics
    path "${output_prefix}.quality_distribution.pdf", emit: quality_distribution_pdf
    path "${output_prefix}.quality_distribution_metrics", emit: quality_distribution_metrics

    script:
    """
    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
        CollectMultipleMetrics \
        --INPUT ${bam} \
        --OUTPUT ${output_prefix} \
        --ASSUME_SORTED true \
        --PROGRAM null \
        --PROGRAM CollectBaseDistributionByCycle \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM QualityScoreDistribution \
        --METRIC_ACCUMULATION_LEVEL null \
        --METRIC_ACCUMULATION_LEVEL ALL_READS
    """
}

process CROSS_CHECK_FINGERPRINTS {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 3500.MB

    input:
    path bams
    path bam_indices
    path haplotype_database
    val lod_threshold
    val cross_check_by
    val output_prefix

    output:
    path "${output_prefix}.cross_check_fingerprints_metrics.txt", emit: cross_check_fingerprints_metrics

    script:
    def inputs_arg = bams.collect { bam -> "--INPUT ${bam}" }.join(' ')
    """
    java -Dsamjdk.buffer_size=131072 \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx3000m \
        -jar /usr/picard/picard.jar \
        CrosscheckFingerprints \
        --OUTPUT ${output_prefix}.cross_check_fingerprints_metrics.txt \
        --HAPLOTYPE_MAP ${haplotype_database} \
        --EXPECT_ALL_GROUPS_TO_MATCH true \
        ${inputs_arg} \
        --LOD_THRESHOLD ${lod_threshold} \
        --CROSSCHECK_BY ${cross_check_by}
    """
}

process GATHER_BAM_FILES {
    
    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 3000.MB

    input:
    path bams
    val create_index
    val compression_level
    val output_prefix

    output:
    path "${output_prefix}.aggregated.bam", emit: aggregated_bam
    path "${output_prefix}.aggregated.bai", emit: aggregated_bam_index, optional: true
    path "${output_prefix}.aggregated.bam.md5", emit: aggregated_bam_md5, optional: true

    script:
    def java_initial_memory_mb = task.memory.toMega() - 1000
    def java_max_memory_mb = task.memory.toMega() - 500
    """
    # Order files numerically using the indices in the filename to get a sorted output.
    bash_inputs_arg=\$(ls ${bams.join(' ')} | sort -V | awk '{print "--INPUT "\$0}' | tr '\\n' ' ')

    java -Dsamjdk.compression_level=${compression_level} \
        -Xms${java_initial_memory_mb}m -Xmx${java_max_memory_mb}m -jar /usr/picard/picard.jar \
        GatherBamFiles \
        \$bash_inputs_arg \
        --OUTPUT ${output_prefix}.aggregated.bam \
        --CREATE_INDEX ${create_index} \
        --CREATE_MD5_FILE ${create_index}
    """
}

process MARK_DUPLICATES {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 7.GB

    input:
    path bams
    val read_name_regex
    val sorting_collection_size_ratio
    val compression_level
    val output_prefix

    output:
    path "${output_prefix}.duplicates_marked.bam", emit: dup_marked_bam
    path "${output_prefix}.duplicate_metrics.txt", emit: duplicate_metrics

    script:
    def inputs_arg = bams.collect { bam -> "--INPUT ${bam}" }.join(' ')
    def read_name_regex_arg = read_name_regex ? "--READ_NAME_REGEX ${read_name_regex}" : ""
    def sorting_collection_size_ratio_arg = sorting_collection_size_ratio ? "--SORTING_COLLECTION_SIZE_RATIO ${sorting_collection_size_ratio}" : ""
    def java_memory_size = task.memory.toGiga() - 2
    """
    java -Dsamjdk.compression_level=${compression_level} -Xms${java_memory_size}g -jar /usr/picard/picard.jar \
        MarkDuplicates \
        ${inputs_arg} \
        --OUTPUT ${output_prefix}.duplicates_marked.bam \
        --METRICS_FILE ${output_prefix}.duplicate_metrics.txt \
        --VALIDATION_STRINGENCY SILENT \
        ${read_name_regex_arg} \
        ${sorting_collection_size_ratio_arg} \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER "queryname" \
        --CLEAR_DT "false" \
        --ADD_PG_TAG_TO_READS false
    """
}

process SORT_SAM {

    container 'us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10'
    memory 5000.MB

    input:
    path bam
    val compression_level
    val output_prefix

    output:
    path "${output_prefix}.sorted.bam", emit: sorted_bam
    path "${output_prefix}.sorted.bai", emit: sorted_bam_index
    path "${output_prefix}.sorted.bam.md5", emit: sorted_bam_md5

    script:
    def java_inital_memory_mb = task.memory.toMega() - 1000
    def java_max_memory_mb = task.memory.toMega() - 500
    """
    java -Dsamjdk.compression_level=${compression_level} \
        -Xms${java_inital_memory_mb}m -Xmx${java_max_memory_mb}m -jar /usr/picard/picard.jar \
        SortSam \
        --INPUT ${bam} \
        --OUTPUT ${output_prefix}.sorted.bam \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE true \
        --MAX_RECORDS_IN_RAM 300000
    """
}
