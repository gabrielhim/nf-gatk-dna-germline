#!/usr/bin/env nextflow

include { AGGREGATE_BAM_QC } from '../../subworkflows/aggregate_bam_qc'
include { BAM_TO_CRAM } from '../../subworkflows/bam_to_cram'
include { UNMAPPED_BAM_TO_ALIGNED_BAM } from '../../subworkflows/unmapped_bam_to_aligned_bam'
include { VARIANT_CALLING } from '../../subworkflows/variant_calling'

include { COLLECT_RAW_WGS_METRICS; COLLECT_WGS_METRICS } from '../../modules/picard'

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

    reference_bin_file = file(params.reference_bin)
    hash_table_cfg_bin_file = file(params.hash_table_cfg_bin)
    hash_table_cmp_file = file(params.hash_table_cmp)
    ref_str_file = file(params.ref_str)

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
    wgs_coverage_interval_list_file = file(params.wgs_coverage_interval_list)

    if (params.dragen_functional_equivalence_mode && params.dragen_maximum_quality_mode) {
        error(
            """
            Both dragen_functional_equivalence_mode and dragen_maximum_quality_mode have been set to true, however, they are mutually exclusive.
            You can set either of them to true, or set them both to false and adjust the arguments individually.
            """.stripIndent()
        )
    }

    run_dragen_mode_variant_calling = params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode ? true : params.run_dragen_mode_variant_calling
    use_spanning_event_genotyping = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.use_spanning_event_genotyping)
    unmap_contaminant_reads = params.dragen_functional_equivalence_mode ? false : (params.dragen_maximum_quality_mode ? true : params.unmap_contaminant_reads)
    perform_bqsr = params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode ? false : params.perform_bqsr
    use_bwa_mem = params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode ? false : params.use_bwa_mem
    use_dragen_hard_filtering = params.dragen_functional_equivalence_mode || params.dragen_maximum_quality_mode ? true : params.use_dragen_hard_filtering

    lod_threshold = -20.0
    cross_check_fingerprints_by = "READGROUP"
    collect_gc_bias_metrics = true
    read_length = 250

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
        reference_bin_file,
        hash_table_cfg_bin_file,
        hash_table_cmp_file,
        dbsnp_vcf_file,
        dbsnp_vcf_index_file,
        known_indels_sites_vcf_files,
        known_indels_sites_vcf_index_files,
        contamination_sites_ud_file,
        contamination_sites_bed_file,
        contamination_sites_mu_file,
        haplotype_database_file,
        params.hard_clip_reads,
        unmap_contaminant_reads,
        params.bin_base_qualities,
        perform_bqsr,
        use_bwa_mem,
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

    wgs_metrics_ch = COLLECT_WGS_METRICS(
        alignment_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        wgs_coverage_interval_list_file,
        read_length,
        params.sample_name,
    )

    wgs_raw_metrics_ch = COLLECT_RAW_WGS_METRICS(
        alignment_ch.final_bam,
        ref_fasta_file,
        ref_fasta_index_file,
        wgs_coverage_interval_list_file,
        read_length,
        params.sample_name,
    )

    variant_calling_ch = VARIANT_CALLING(
        params.sample_name,
        alignment_ch.final_bam,
        calling_interval_list_file,
        evaluation_interval_list_file,
        ref_fasta_file,
        ref_fasta_index_file,
        ref_dict_file,
        ref_str_file,
        dbsnp_vcf_file,
        dbsnp_vcf_index_file,
        params.haplotype_scatter_count,
        params.break_bands_at_multiples_of,
        alignment_ch.contamination,
        run_dragen_mode_variant_calling,
        use_spanning_event_genotyping,
        true,
        params.save_bamout,
        params.skip_reblocking,
        use_dragen_hard_filtering,
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

    wgs_metrics = wgs_metrics_ch.wgs_metrics
    wgs_raw_metrics = wgs_raw_metrics_ch.wgs_raw_metrics

    final_vcf = variant_calling_ch.final_vcf.map { vcf, _vcf_index -> vcf }
    final_vcf_index = variant_calling_ch.final_vcf.map { _vcf, vcf_index -> vcf_index }
    bamout = variant_calling_ch.bamout.map { bam, _bam_index -> bam }
    bamout_index = variant_calling_ch.bamout.map { _bam, bam_index -> bam_index }

    vcf_summary_metrics = variant_calling_ch.vcf_summary_metrics
    vcf_detail_metrics = variant_calling_ch.vcf_detail_metrics
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
    wgs_metrics {
        path 'quality_control'
    }
    wgs_raw_metrics {
        path 'quality_control'
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
}
