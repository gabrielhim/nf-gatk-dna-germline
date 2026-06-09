#!/usr/bin/env nextflow

include { AGGREGATE_BAM_QC } from '../../subworkflows/aggregate_bam_qc'
include { BAM_TO_CRAM } from '../../subworkflows/bam_to_cram'
include { UNMAPPED_BAM_TO_ALIGNED_BAM } from '../../subworkflows/unmapped_bam_to_aligned_bam'
include { VARIANT_CALLING } from '../../subworkflows/variant_calling'

include { SUBSET_CONTAMINATION_RESOURCES } from '../../modules/bedtools'
include { COLLECT_HS_METRICS } from '../../modules/picard'

workflow {
    
    main:
    unmapped_bams_ch = channel.fromPath(params.unmapped_bams).map { ubam -> [ubam, ubam.baseName] }

    ref_fasta_file = file(params.ref_fasta)
    ref_fasta_index_file = file(params.ref_fasta_index)
    ref_dict_file = file(params.ref_dict)

    ref_alt_file = file(params.ref_alt)
    ref_amb_file = file(params.ref_amb)
    ref_ann_file = file(params.ref_ann)
    ref_bwt_file = file(params.ref_bwt)
    ref_pac_file = file(params.ref_pac)
    ref_sa_file = file(params.ref_sa)

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

    calling_interval_list_file = file(params.calling_interval_list)
    evaluation_interval_list_file = file(params.evaluation_interval_list)
    target_interval_list_file = file(params.target_interval_list)
    bait_interval_list_file = file(params.bait_interval_list)

    lod_threshold = -10.0
    cross_check_fingerprints_by = "READGROUP"
    collect_gc_bias_metrics = false
    make_gvcf = true

    subset_contamination_sites_ch = SUBSET_CONTAMINATION_RESOURCES(
        target_interval_list_file,
        contamination_sites_ud_file,
        contamination_sites_bed_file,
        contamination_sites_mu_file,
    )

    alignment_ch = UNMAPPED_BAM_TO_ALIGNED_BAM(
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
        [],
        [],
        [],
        dbsnp_vcf_file,
        dbsnp_vcf_index_file,
        known_indels_sites_vcf_files,
        known_indels_sites_vcf_index_files,
        subset_contamination_sites_ch.subset_ud,
        subset_contamination_sites_ch.subset_bed,
        subset_contamination_sites_ch.subset_mu,
        haplotype_database_file,
        params.hard_clip_reads,
        params.unmap_contaminant_reads,
        params.bin_base_qualities,
        params.perform_bqsr,
        params.use_bwa_mem,
        params.allow_empty_ref_alt,
        lod_threshold,
        cross_check_fingerprints_by,
    )

    aggr_bam_qc_ch = AGGREGATE_BAM_QC(
        params.sample_name,
        alignment_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        haplotype_database_file,
        fingerprint_genotypes_file,
        fingerprint_genotypes_index_file,
        collect_gc_bias_metrics,
    )

    bam_to_cram_ch = BAM_TO_CRAM(
        params.sample_name,
        alignment_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        alignment_ch.duplication_metrics,
        aggr_bam_qc_ch.aggr_alignment_summary_metrics,
    )

    variant_calling_ch = VARIANT_CALLING(
        params.sample_name,
        alignment_ch.final_bam,
        calling_interval_list_file,
        evaluation_interval_list_file,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        [],
        dbsnp_vcf_file,
        dbsnp_vcf_index_file,
        params.haplotype_scatter_count,
        params.break_bands_at_multiples_of,
        alignment_ch.contamination,
        params.run_dragen_mode_variant_calling,
        params.use_spanning_event_genotyping,
        make_gvcf,
        params.save_bamout,
        params.skip_reblocking,
        params.use_dragen_hard_filtering,
    )

    hs_metrics_ch = COLLECT_HS_METRICS(
        alignment_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        target_interval_list_file,
        bait_interval_list_file,
        params.sample_name,
    )

    publish:
    quality_yield_metrics = alignment_ch.quality_yield_metrics

    unsorted_read_group_base_distribution_by_cycle_pdf = alignment_ch.unsorted_read_group_base_distribution_by_cycle_pdf
    unsorted_read_group_base_distribution_by_cycle_metrics = alignment_ch.unsorted_read_group_base_distribution_by_cycle_metrics
    unsorted_read_group_insert_size_histogram_pdf = alignment_ch.unsorted_read_group_insert_size_histogram_pdf
    unsorted_read_group_insert_size_metrics = alignment_ch.unsorted_read_group_insert_size_metrics
    unsorted_read_group_quality_by_cycle_pdf = alignment_ch.unsorted_read_group_quality_by_cycle_pdf
    unsorted_read_group_quality_by_cycle_metrics = alignment_ch.unsorted_read_group_quality_by_cycle_metrics
    unsorted_read_group_quality_distribution_pdf = alignment_ch.unsorted_read_group_quality_distribution_pdf
    unsorted_read_group_quality_distribution_metrics = alignment_ch.unsorted_read_group_quality_distribution_metrics

    cross_check_fingerprints_metrics = alignment_ch.cross_check_fingerprints_metrics

    per_sample_stats = alignment_ch.per_sample_stats
    contamination = alignment_ch.contamination

    duplication_metrics = alignment_ch.duplication_metrics
    bqsr_report = alignment_ch.bqsr_report

    final_bam = alignment_ch.final_bam.map { bam, _bai -> bam }
    final_bam_index = alignment_ch.final_bam.map { _bam, bai -> bai }

    read_group_alignment_summary_metrics = aggr_bam_qc_ch.read_group_alignment_summary_metrics
    read_group_gc_bias_detail_metrics = aggr_bam_qc_ch.read_group_gc_bias_detail_metrics
    read_group_gc_bias_pdf = aggr_bam_qc_ch.read_group_gc_bias_pdf
    read_group_gc_bias_summary_metrics = aggr_bam_qc_ch.read_group_gc_bias_summary_metrics

    read_group_checksum = aggr_bam_qc_ch.read_group_checksum

    aggr_alignment_summary_metrics = aggr_bam_qc_ch.aggr_alignment_summary_metrics
    aggr_bait_bias_detail_metrics = aggr_bam_qc_ch.aggr_bait_bias_detail_metrics
    aggr_bait_bias_summary_metrics = aggr_bam_qc_ch.aggr_bait_bias_summary_metrics
    aggr_gc_bias_detail_metrics = aggr_bam_qc_ch.aggr_gc_bias_detail_metrics
    aggr_gc_bias_pdf = aggr_bam_qc_ch.aggr_gc_bias_pdf
    aggr_gc_bias_summary_metrics = aggr_bam_qc_ch.aggr_gc_bias_summary_metrics
    aggr_insert_size_histogram_pdf = aggr_bam_qc_ch.aggr_insert_size_histogram_pdf
    aggr_insert_size_metrics = aggr_bam_qc_ch.aggr_insert_size_metrics
    aggr_pre_adapter_detail_metrics = aggr_bam_qc_ch.aggr_pre_adapter_detail_metrics
    aggr_pre_adapter_summary_metrics = aggr_bam_qc_ch.aggr_pre_adapter_summary_metrics
    aggr_quality_distribution_pdf = aggr_bam_qc_ch.aggr_quality_distribution_pdf
    aggr_quality_distribution_metrics = aggr_bam_qc_ch.aggr_quality_distribution_metrics
    aggr_error_summary_metrics = aggr_bam_qc_ch.aggr_error_summary_metrics

    fingerprint_summary_metrics = aggr_bam_qc_ch.fingerprint_summary_metrics
    fingerprint_detail_metrics = aggr_bam_qc_ch.fingerprint_detail_metrics

    cram = bam_to_cram_ch.cram.map { cram, _cram_index -> cram }
    cram_index = bam_to_cram_ch.cram.map { _cram, cram_index -> cram_index }
    cram_md5 = bam_to_cram_ch.cram_md5
    cram_validation_report = bam_to_cram_ch.validation_report

    final_vcf = variant_calling_ch.final_vcf.map { vcf, _vcf_index -> vcf }
    final_vcf_index = variant_calling_ch.final_vcf.map { _vcf, vcf_index -> vcf_index }
    bamout = variant_calling_ch.bamout.map { bam, _bam_index -> bam }
    bamout_index = variant_calling_ch.bamout.map { _bam, bam_index -> bam_index }

    vcf_summary_metrics = variant_calling_ch.vcf_summary_metrics
    vcf_detail_metrics = variant_calling_ch.vcf_detail_metrics

    hybrid_selection_metrics = hs_metrics_ch.hybrid_selection_metrics
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
    aggr_alignment_summary_metrics {
        path 'quality_control'
    }
    aggr_bait_bias_detail_metrics {
        path 'quality_control'
    }
    aggr_bait_bias_summary_metrics {
        path 'quality_control'
    }
    aggr_gc_bias_detail_metrics {
        path 'quality_control'
    }
    aggr_gc_bias_pdf {
        path 'quality_control'
    }
    aggr_gc_bias_summary_metrics {
        path 'quality_control'
    }
    aggr_insert_size_histogram_pdf {
        path 'quality_control'
    }
    aggr_insert_size_metrics {
        path 'quality_control'
    }
    aggr_pre_adapter_detail_metrics {
        path 'quality_control'
    }
    aggr_pre_adapter_summary_metrics {
        path 'quality_control'
    }
    aggr_quality_distribution_pdf {
        path 'quality_control'
    }
    aggr_quality_distribution_metrics {
        path 'quality_control'
    }
    aggr_error_summary_metrics {
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
    final_vcf {
        path 'variants'
    }
    final_vcf_index {
        path 'variants'
    }
    bamout {
        path 'variants'
    }
    bamout_index {
        path 'variants'
    }
    vcf_summary_metrics {
        path 'quality_control'
    }
    vcf_detail_metrics {
        path 'quality_control'
    }
    hybrid_selection_metrics {
        path 'quality_control'
    }
}
