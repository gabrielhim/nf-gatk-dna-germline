#!/usr/bin/env nextflow

include { 
    CALCULATE_READ_GROUP_CHECKSUM;
    CHECK_FINGERPRINT_TASK;
    COLLECT_AGGREGATION_METRICS;
    COLLECT_READ_GROUP_BAM_METRICS;
} from '../modules/picard'

workflow AGGREGATE_BAM_QC {

    take:
    sample_name
    bam_ch
    ref_fasta
    ref_fasta_index
    ref_dict
    haplotype_database
    fingerprint_genotypes
    fingerprint_genotypes_index
    collect_gc_bias_metrics

    main:
    genotype_lod_threshold = 5.0

    COLLECT_READ_GROUP_BAM_METRICS(
        bam_ch,
        ref_fasta,
        ref_fasta_index,
        ref_dict,
        collect_gc_bias_metrics,
        "${sample_name}.read_group"
    )

    COLLECT_AGGREGATION_METRICS(
        bam_ch, ref_fasta, ref_fasta_index, ref_dict, collect_gc_bias_metrics, sample_name
    )

    check_fingerprint = haplotype_database && fingerprint_genotypes
    if (check_fingerprint) {
        CHECK_FINGERPRINT_TASK(
            bam_ch,
            channel.empty(),
            fingerprint_genotypes,
            fingerprint_genotypes_index,
            haplotype_database,
            ref_fasta,
            ref_fasta_index,
            genotype_lod_threshold,
            false,
            [],
            sample_name,
            sample_name,
        )
    }

    CALCULATE_READ_GROUP_CHECKSUM(bam_ch, sample_name)

    emit:
    read_group_alignment_summary_metrics = COLLECT_READ_GROUP_BAM_METRICS.out.alignment_summary_metrics
    read_group_gc_bias_detail_metrics = COLLECT_READ_GROUP_BAM_METRICS.out.gc_bias_detail_metrics
    read_group_gc_bias_pdf = COLLECT_READ_GROUP_BAM_METRICS.out.gc_bias_pdf
    read_group_gc_bias_summary_metrics = COLLECT_READ_GROUP_BAM_METRICS.out.gc_bias_summary_metrics

    read_group_checksum = CALCULATE_READ_GROUP_CHECKSUM.out.read_group_md5

    aggr_alignment_summary_metrics = COLLECT_AGGREGATION_METRICS.out.alignment_summary_metrics
    aggr_bait_bias_detail_metrics = COLLECT_AGGREGATION_METRICS.out.bait_bias_detail_metrics
    aggr_bait_bias_summary_metrics = COLLECT_AGGREGATION_METRICS.out.bait_bias_summary_metrics
    aggr_gc_bias_detail_metrics = collect_gc_bias_metrics ? COLLECT_AGGREGATION_METRICS.out.gc_bias_detail_metrics : []
    aggr_gc_bias_pdf = collect_gc_bias_metrics ? COLLECT_AGGREGATION_METRICS.out.gc_bias_pdf : []
    aggr_gc_bias_summary_metrics = collect_gc_bias_metrics ? COLLECT_AGGREGATION_METRICS.out.gc_bias_summary_metrics : []
    aggr_insert_size_histogram_pdf = collect_gc_bias_metrics ? COLLECT_AGGREGATION_METRICS.out.insert_size_histogram_pdf : []
    aggr_insert_size_metrics = collect_gc_bias_metrics ? COLLECT_AGGREGATION_METRICS.out.insert_size_metrics : []
    aggr_pre_adapter_detail_metrics = COLLECT_AGGREGATION_METRICS.out.pre_adapter_detail_metrics
    aggr_pre_adapter_summary_metrics = COLLECT_AGGREGATION_METRICS.out.pre_adapter_summary_metrics
    aggr_quality_distribution_pdf = COLLECT_AGGREGATION_METRICS.out.quality_distribution_pdf
    aggr_quality_distribution_metrics = COLLECT_AGGREGATION_METRICS.out.quality_distribution_metrics
    aggr_error_summary_metrics = COLLECT_AGGREGATION_METRICS.out.error_summary_metrics

    fingerprint_summary_metrics = check_fingerprint ? CHECK_FINGERPRINT_TASK.out.summary_metrics : []
    fingerprint_detail_metrics = check_fingerprint ? CHECK_FINGERPRINT_TASK.out.detail_metrics : []
}
