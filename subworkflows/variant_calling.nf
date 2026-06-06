#!/usr/bin/env nextflow

include { 
    CALIBRATE_DRAGSTR_MODEL;
    DRAGEN_HARD_FILTER_VCF;
    HAPLOTYPE_CALLER;
    REBLOCK_GVCF;
    VALIDATE_VCF;
} from '../modules/gatk'
include {
    COLLECT_VARIANT_CALLING_METRICS;
    MERGE_VCFS;
    SCATTER_INTERVAL_LIST;
    SORT_SAM;
} from '../modules/picard'
include { MERGE_BAMOUTS } from '../modules/samtools'

workflow VARIANT_CALLING {
    
    take:
    sample_name
    bam_ch
    bam_index_ch
    calling_interval_list
    evaluation_interval_list
    ref_fasta
    ref_fasta_index
    ref_dict
    ref_str

    dbsnp_vcf
    dbsnp_vcf_index

    haplotype_scatter_count
    break_bands_at_multiples_of
    contamination_ch

    run_dragen_mode_variant_calling
    use_spanning_event_genotyping
    make_gvcf
    make_bamout
    skip_reblocking
    use_dragen_hard_filtering

    main:
    if (run_dragen_mode_variant_calling) {
        CALIBRATE_DRAGSTR_MODEL(
            bam_ch, bam_index_ch, ref_fasta, ref_fasta_index, ref_dict, ref_str, sample_name
        )
    }

    SCATTER_INTERVAL_LIST(
        calling_interval_list, haplotype_scatter_count, break_bands_at_multiples_of
    )

    haplotype_caller_ch = HAPLOTYPE_CALLER(
        bam_ch,
        bam_index_ch,
        SCATTER_INTERVAL_LIST.out.scattered_interval_lists.flatten(),
        ref_fasta,
        ref_fasta_index,
        ref_dict,
        run_dragen_mode_variant_calling ? CALIBRATE_DRAGSTR_MODEL.out.dragstr_model : [],
        make_gvcf,
        make_bamout,
        run_dragen_mode_variant_calling,
        use_spanning_event_genotyping,
        run_dragen_mode_variant_calling ? 0 : contamination_ch,
        sample_name,
    )

    merged_vcf_ch = MERGE_VCFS(
        haplotype_caller_ch.vcf.collect(),
        haplotype_caller_ch.vcf_index.collect(),
        make_gvcf ? "${sample_name}.haplotype_caller.g.vcf.gz" : "${sample_name}.haplotype_caller.vcf.gz",
    )

    if (make_bamout) {
        SORT_SAM(haplotype_caller_ch.bamout, 2, sample_name, true)

        MERGE_BAMOUTS(SORT_SAM.out.sorted_bam.collect(), "${sample_name}.haplotype_caller.bamout.bam")
    }

    if (use_dragen_hard_filtering) {
        DRAGEN_HARD_FILTER_VCF(
            merged_vcf_ch.merged_vcf,
            merged_vcf_ch.merged_vcf_index,
            make_gvcf ? "${sample_name}.hard-filtered.g.vcf.gz" : "${sample_name}.hard-filtered.vcf.gz",
        )
    }

    merged_filt_vcf_ch = use_dragen_hard_filtering ? DRAGEN_HARD_FILTER_VCF.out.filtered_vcf : merged_vcf_ch.merged_vcf
    merged_filt_vcf_index_ch = use_dragen_hard_filtering ? DRAGEN_HARD_FILTER_VCF.out.filtered_vcf_index : merged_vcf_ch.merged_vcf_index

    run_reblock = make_gvcf && !skip_reblocking
    if (run_reblock) {
        REBLOCK_GVCF(
            merged_filt_vcf_ch,
            merged_filt_vcf_index_ch,
            ref_fasta,
            ref_fasta_index,
            ref_dict,
            [],
            false,
            sample_name,
        )
    }

    final_vcf = run_reblock ? REBLOCK_GVCF.out.reblocked_gvcf : merged_filt_vcf_ch
    final_vcf_index = run_reblock ? REBLOCK_GVCF.out.reblocked_gvcf_index : merged_filt_vcf_index_ch

    VALIDATE_VCF(
        final_vcf,
        final_vcf_index,
        ref_fasta,
        ref_fasta_index,
        ref_dict,
        dbsnp_vcf,
        dbsnp_vcf_index,
        calling_interval_list,
        make_gvcf,
        !skip_reblocking ? "--no-overlaps" : "",
    )

    COLLECT_VARIANT_CALLING_METRICS(
        final_vcf,
        final_vcf_index,
        ref_dict,
        dbsnp_vcf,
        dbsnp_vcf_index,
        evaluation_interval_list,
        make_gvcf,
        sample_name,
    )

    emit:
    final_vcf = final_vcf
    final_vcf_index = final_vcf_index
    bamout = make_bamout ? MERGE_BAMOUTS.out.merged_bamout : []
    bamout_index = make_bamout ? MERGE_BAMOUTS.out.merged_bamout_index : []
    vcf_summary_metrics = COLLECT_VARIANT_CALLING_METRICS.out.variant_calling_summary_metrics
    vcf_detail_metrics = COLLECT_VARIANT_CALLING_METRICS.out.variant_calling_detail_metrics
}
