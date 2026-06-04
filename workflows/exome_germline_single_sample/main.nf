#!/usr/bin/env nextflow

include { SUBSET_CONTAMINATION_RESOURCES } from '../../modules/bedtools'

include { AGGREGATE_BAM_QC } from '../../subworkflows/aggregate_bam_qc'
include { BAM_TO_CRAM } from '../../subworkflows/bam_to_cram'
include { UNMAPPED_BAM_TO_ALIGNED_BAM } from '../../subworkflows/unmapped_bam_to_aligned_bam'

workflow {
    
    main:
    intervals_ch = channel.fromPath(params.target_interval_list)

    ref_fasta_file = file(params.ref_fasta)
    ref_fasta_index_file = file(params.ref_fasta_index)
    ref_dict_file = file(params.ref_dict)

    ref_alt_file = file(params.ref_alt)
    ref_amb_file = file(params.ref_amb)
    ref_ann_file = file(params.ref_ann)
    ref_bwt_file = file(params.ref_bwt)
    ref_pac_file = file(params.ref_pac)
    ref_sa_file = file(params.ref_sa)

    reference_bin_file = file(params.reference_bin)
    hash_table_cfg_bin_file = file(params.hash_table_cfg_bin)
    hash_table_cmp_file = file(params.hash_table_cmp)

    dbsnp_vcf_file = file(params.dbsnp_vcf)
    dbsnp_vcf_index_file = file(params.dbsnp_vcf_index)
    known_indels_sites_vcf_files = channel.fromPath(params.known_indels_sites_vcfs).collect()
    known_indels_sites_vcf_index_files = channel.fromPath(params.known_indels_sites_indices).collect()

    contamination_sites_ud_file = file(params.contamination_sites_ud)
    contamination_sites_bed_file = file(params.contamination_sites_bed)
    contamination_sites_mu_file = file(params.contamination_sites_mu)
    haplotype_database_file = file(params.haplotype_database)
    fingerprint_genotypes_file = file(params.fingerprint_genotypes)
    fingerprint_genotypes_index_file = file(params.fingerprint_genotypes_index)

    subset_contamination_sites_ch = SUBSET_CONTAMINATION_RESOURCES(
        intervals_ch,
        contamination_sites_ud_file,
        contamination_sites_bed_file,
        contamination_sites_mu_file,
    )

    unmapped_bams_ch = channel.fromPath(params.unmapped_bams).map { ubam -> [ubam, ubam.baseName] }

    mapping_workflow_ch = UNMAPPED_BAM_TO_ALIGNED_BAM(
        params.sample_name,
        unmapped_bams_ch,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        ref_alt_file,
        ref_amb_file,
        ref_ann_file,
        ref_bwt_file,
        ref_pac_file,
        ref_sa_file,
        reference_bin_file,
        hash_table_cfg_bin_file,
        hash_table_cmp_file,
        dbsnp_vcf_file,
        dbsnp_vcf_index_file,
        known_indels_sites_vcf_files,
        known_indels_sites_vcf_index_files,
        subset_contamination_sites_ch.subset_ud,
        subset_contamination_sites_ch.subset_bed,
        subset_contamination_sites_ch.subset_mu,
        haplotype_database_file,
        params.use_bwa_mem,
        params.perform_bqsr,
    )

    agg_bam_qc_ch = AGGREGATE_BAM_QC(
        params.sample_name,
        mapping_workflow_ch.final_bam,
        mapping_workflow_ch.final_bam_index,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        haplotype_database_file,
        fingerprint_genotypes_file,
        fingerprint_genotypes_index_file,
    )

    bam_to_cram_ch = BAM_TO_CRAM(
        params.sample_name,
        mapping_workflow_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        mapping_workflow_ch.duplication_metrics,
        agg_bam_qc_ch.agg_alignment_summary_metrics,
    )

    publish:
    quality_yield_metrics = mapping_workflow_ch.quality_yield_metrics

    unsorted_read_group_base_distribution_by_cycle_pdf = mapping_workflow_ch.unsorted_read_group_base_distribution_by_cycle_pdf
    unsorted_read_group_base_distribution_by_cycle_metrics = mapping_workflow_ch.unsorted_read_group_base_distribution_by_cycle_metrics
    unsorted_read_group_insert_size_histogram_pdf = mapping_workflow_ch.unsorted_read_group_insert_size_histogram_pdf
    unsorted_read_group_insert_size_metrics = mapping_workflow_ch.unsorted_read_group_insert_size_metrics
    unsorted_read_group_quality_by_cycle_pdf = mapping_workflow_ch.unsorted_read_group_quality_by_cycle_pdf
    unsorted_read_group_quality_by_cycle_metrics = mapping_workflow_ch.unsorted_read_group_quality_by_cycle_metrics
    unsorted_read_group_quality_distribution_pdf = mapping_workflow_ch.unsorted_read_group_quality_distribution_pdf
    unsorted_read_group_quality_distribution_metrics = mapping_workflow_ch.unsorted_read_group_quality_distribution_metrics

    cross_check_fingerprints_metrics = mapping_workflow_ch.cross_check_fingerprints_metrics

    per_sample_stats = mapping_workflow_ch.per_sample_stats
    contamination = mapping_workflow_ch.contamination

    duplication_metrics = mapping_workflow_ch.duplication_metrics
    bqsr_report = mapping_workflow_ch.bqsr_report

    final_bam = mapping_workflow_ch.final_bam
    final_bam_index = mapping_workflow_ch.final_bam_index

    read_group_alignment_summary_metrics = agg_bam_qc_ch.read_group_alignment_summary_metrics
    read_group_gc_bias_detail_metrics = agg_bam_qc_ch.read_group_gc_bias_detail_metrics
    read_group_gc_bias_pdf = agg_bam_qc_ch.read_group_gc_bias_pdf
    read_group_gc_bias_summary_metrics = agg_bam_qc_ch.read_group_gc_bias_summary_metrics

    read_group_checksum = agg_bam_qc_ch.read_group_checksum

    agg_alignment_summary_metrics = agg_bam_qc_ch.agg_alignment_summary_metrics
    agg_bait_bias_detail_metrics = agg_bam_qc_ch.agg_bait_bias_detail_metrics
    agg_bait_bias_summary_metrics = agg_bam_qc_ch.agg_bait_bias_summary_metrics
    agg_gc_bias_detail_metrics = agg_bam_qc_ch.agg_gc_bias_detail_metrics
    agg_gc_bias_pdf = agg_bam_qc_ch.agg_gc_bias_pdf
    agg_gc_bias_summary_metrics = agg_bam_qc_ch.agg_gc_bias_summary_metrics
    agg_insert_size_histogram_pdf = agg_bam_qc_ch.agg_insert_size_histogram_pdf
    agg_insert_size_metrics = agg_bam_qc_ch.agg_insert_size_metrics
    agg_pre_adapter_detail_metrics = agg_bam_qc_ch.agg_pre_adapter_detail_metrics
    agg_pre_adapter_summary_metrics = agg_bam_qc_ch.agg_pre_adapter_summary_metrics
    agg_quality_distribution_pdf = agg_bam_qc_ch.agg_quality_distribution_pdf
    agg_quality_distribution_metrics = agg_bam_qc_ch.agg_quality_distribution_metrics
    agg_error_summary_metrics = agg_bam_qc_ch.agg_error_summary_metrics

    fingerprint_summary_metrics = agg_bam_qc_ch.fingerprint_summary_metrics
    fingerprint_detail_metrics = agg_bam_qc_ch.fingerprint_detail_metrics

    cram = bam_to_cram_ch.cram
    cram_index = bam_to_cram_ch.cram_index
    cram_md5 = bam_to_cram_ch.cram_md5
    cram_validation_report = bam_to_cram_ch.validation_report
}

output {
    quality_yield_metrics {
        path 'quality_control'
    }
    unsorted_read_group_base_distribution_by_cycle_pdf {
        path 'quality_control'
    }
    unsorted_read_group_base_distribution_by_cycle_metrics {
        path 'quality_control'
    }
    unsorted_read_group_insert_size_histogram_pdf {
        path 'quality_control'
    }
    unsorted_read_group_insert_size_metrics {
        path 'quality_control'
    }
    unsorted_read_group_quality_by_cycle_pdf {
        path 'quality_control'
    }
    unsorted_read_group_quality_by_cycle_metrics {
        path 'quality_control'
    }
    unsorted_read_group_quality_distribution_pdf {
        path 'quality_control'
    }
    unsorted_read_group_quality_distribution_metrics {
        path 'quality_control'
    }
    cross_check_fingerprints_metrics {
        path 'quality_control'
    }
    per_sample_stats {
        path 'quality_control'
    }
    contamination {
        path 'quality_control'
    }
    duplication_metrics {
        path 'duplicates'
    }
    bqsr_report {
        path 'recalibration'
    }
    final_bam {
        path 'alignment'
    }
    final_bam_index {
        path 'alignment'
    }
    read_group_alignment_summary_metrics {
        path 'quality_control'
    }
    read_group_gc_bias_detail_metrics {
        path 'quality_control'
    }
    read_group_gc_bias_pdf {
        path 'quality_control'
    }
    read_group_gc_bias_summary_metrics {
        path 'quality_control'
    }
    read_group_checksum {
        path 'quality_control'
    }
    agg_alignment_summary_metrics {
        path 'quality_control'
    }
    agg_bait_bias_detail_metrics {
        path 'quality_control'
    }
    agg_bait_bias_summary_metrics {
        path 'quality_control'
    }
    agg_gc_bias_detail_metrics {
        path 'quality_control'
    }
    agg_gc_bias_pdf {
        path 'quality_control'
    }
    agg_gc_bias_summary_metrics {
        path 'quality_control'
    }
    agg_insert_size_histogram_pdf {
        path 'quality_control'
    }
    agg_insert_size_metrics {
        path 'quality_control'
    }
    agg_pre_adapter_detail_metrics {
        path 'quality_control'
    }
    agg_pre_adapter_summary_metrics {
        path 'quality_control'
    }
    agg_quality_distribution_pdf {
        path 'quality_control'
    }
    agg_quality_distribution_metrics {
        path 'quality_control'
    }
    agg_error_summary_metrics {
        path 'quality_control'
    }
    fingerprint_summary_metrics {
        path 'quality_control'
    }
    fingerprint_detail_metrics {
        path 'quality_control'
    }
    cram {
        path 'alignment'
    }
    cram_index {
        path 'alignment'
    }
    cram_md5 {
        path 'alignment'
    }
    cram_validation_report {
        path 'alignment'
    }
}
