#!/usr/bin/env nextflow

include { BWA_MEM_AND_MERGE_ALIGNMENT } from '../modules/bwa'
include { DRAGMAP_AND_MERGE_ALIGNMENT } from '../modules/dragmap'
include { CREATE_SEQUENCE_GROUPING_TSV } from '../modules/python'
include { CHECK_CONTAMINATION } from '../modules/verify_bam_id'
include { APPLY_BQSR; BASE_RECALIBRATOR; GATHER_BQSR_REPORTS } from '../modules/gatk'
include { 
    COLLECT_QUALITY_YIELD_METRICS;
    COLLECT_READ_GROUP_BAM_METRICS;
    COLLECT_UNSORTED_READ_GROUP_BAM_METRICS;
    CROSS_CHECK_FINGERPRINTS;
    GATHER_BAM_FILES;
    MARK_DUPLICATES;
    SORT_SAM;
} from '../modules/picard'
 
workflow UNMAPPED_BAM_TO_ALIGNED_BAM {

    take:
    sample_name
    unmapped_bams_ch
    ref_fasta
    ref_fasta_index
    ref_dict
    ref_alt
    ref_amb
    ref_ann
    ref_bwt
    ref_pac
    ref_sa

    reference_bin
    hash_table_cfg_bin
    hash_table_cmp

    dbsnp_vcf
    dbsnp_vcf_index
    known_indels_sites_vcfs
    known_indels_sites_vcf_indices

    target_contamination_sites_ud
    target_contamination_sites_bed
    target_contamination_sites_mu
    
    haplotype_database

    use_bwa_mem
    
    main:
    compression_level = 2
    
    hard_clip_reads = false
    unmap_contaminant_reads = true
    bin_base_qualities = true
    somatic = false
    perform_bqsr = true
    allow_empty_ref_alt = false

    lod_threshold = -10.0
    cross_check_fingerprints_by = "READGROUP"
    contamination_underestimation_factor = 0.75

    COLLECT_QUALITY_YIELD_METRICS(unmapped_bams_ch)

    // Scatter by unmapped BAM.
    if (use_bwa_mem) {
        BWA_MEM_AND_MERGE_ALIGNMENT(
            unmapped_bams_ch,
            ref_fasta,
            ref_fasta_index,
            ref_dict,
            ref_alt,
            ref_amb,
            ref_ann,
            ref_bwt,
            ref_pac,
            ref_sa,
            compression_level,
            allow_empty_ref_alt,
            hard_clip_reads,
            unmap_contaminant_reads,
        )

        merged_bam_ch = BWA_MEM_AND_MERGE_ALIGNMENT.out.merged_bam
        mapping_stderr_log_ch = BWA_MEM_AND_MERGE_ALIGNMENT.out.bwa_stderr_log

    } else {
        DRAGMAP_AND_MERGE_ALIGNMENT(
            unmapped_bams_ch,
            ref_fasta,
            ref_fasta_index,
            ref_dict,
            reference_bin,
            hash_table_cfg_bin,
            hash_table_cmp,
            compression_level,
            hard_clip_reads,
            unmap_contaminant_reads,
        )

        merged_bam_ch = DRAGMAP_AND_MERGE_ALIGNMENT.out.merged_bam
        mapping_stderr_log_ch = DRAGMAP_AND_MERGE_ALIGNMENT.out.dragmap_stderr_log
    }

    read_group_metrics_input_ch = merged_bam_ch.map { bam -> [ "${bam.baseName}.read_group", bam] }
    COLLECT_UNSORTED_READ_GROUP_BAM_METRICS(read_group_metrics_input_ch)

    MARK_DUPLICATES(merged_bam_ch.collect(), false, false, compression_level, sample_name)

    sorted_bam_ch = SORT_SAM(MARK_DUPLICATES.out.dup_marked_bam, compression_level, sample_name)

    if (haplotype_database) {
        CROSS_CHECK_FINGERPRINTS(
            sorted_bam_ch.sorted_bam,
            sorted_bam_ch.sorted_bam_index,
            haplotype_database,
            lod_threshold,
            cross_check_fingerprints_by,
            sample_name,
        )
    }

    CHECK_CONTAMINATION(
        sorted_bam_ch.sorted_bam,
        sorted_bam_ch.sorted_bam_index,
        target_contamination_sites_ud,
        target_contamination_sites_bed,
        target_contamination_sites_mu,
        ref_fasta,
        ref_fasta_index,
        contamination_underestimation_factor,
        compression_level,
        true,
        sample_name,
    )

    if (perform_bqsr) {
        CREATE_SEQUENCE_GROUPING_TSV(ref_dict)

        sequence_grouping_ch = CREATE_SEQUENCE_GROUPING_TSV.out.sequence_grouping
            .splitText()
            .map { line -> line.trim().split('\t').toList() }
        
        sequence_grouping_with_unmapped_ch = CREATE_SEQUENCE_GROUPING_TSV.out.sequence_grouping_with_unmapped
            .splitText()
            .map { line -> line.trim().split('\t').toList() }

        BASE_RECALIBRATOR(
            sorted_bam_ch.sorted_bam,
            sorted_bam_ch.sorted_bam_index,
            dbsnp_vcf,
            dbsnp_vcf_index,
            known_indels_sites_vcfs,
            known_indels_sites_vcf_indices,
            ref_fasta,
            ref_fasta_index,
            ref_dict,
            sequence_grouping_ch,
            sample_name,
        )

        GATHER_BQSR_REPORTS(BASE_RECALIBRATOR.out.recalibration_report.collect(), sample_name)

        APPLY_BQSR(
            sorted_bam_ch.sorted_bam,
            sorted_bam_ch.sorted_bam_index,
            GATHER_BQSR_REPORTS.out.bqsr_report,
            ref_fasta,
            ref_fasta_index,
            ref_dict,
            sequence_grouping_with_unmapped_ch,
            compression_level,
            bin_base_qualities,
            somatic,
            sample_name,
        )

        GATHER_BAM_FILES(APPLY_BQSR.out.recalibrated_bam.collect(), true, compression_level, sample_name)
    }

    emit:
    quality_yield_metrics = COLLECT_QUALITY_YIELD_METRICS.out.quality_yield_metrics.collect()

    unsorted_read_group_base_distribution_by_cycle_pdf = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.base_distribution_by_cycle_pdf
    unsorted_read_group_base_distribution_by_cycle_metrics = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.base_distribution_by_cycle_metrics
    unsorted_read_group_insert_size_histogram_pdf = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.insert_size_histogram_pdf
    unsorted_read_group_insert_size_metrics = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.insert_size_metrics
    unsorted_read_group_quality_by_cycle_pdf = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.quality_by_cycle_pdf
    unsorted_read_group_quality_by_cycle_metrics = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.quality_by_cycle_metrics
    unsorted_read_group_quality_distribution_pdf = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.quality_distribution_pdf
    unsorted_read_group_quality_distribution_metrics = COLLECT_UNSORTED_READ_GROUP_BAM_METRICS.out.quality_distribution_metrics

    cross_check_fingerprints_metrics = CROSS_CHECK_FINGERPRINTS.out.cross_check_fingerprints_metrics

    per_sample_stats = CHECK_CONTAMINATION.out.per_sample_stats
    contamination = CHECK_CONTAMINATION.out.contamination

    duplicate_metrics = MARK_DUPLICATES.out.duplicate_metrics
    bqsr_report = GATHER_BQSR_REPORTS.out.bqsr_report

    final_bam = GATHER_BAM_FILES.out.aggregated_bam
    final_bam_index = GATHER_BAM_FILES.out.aggregated_bam_index
}
