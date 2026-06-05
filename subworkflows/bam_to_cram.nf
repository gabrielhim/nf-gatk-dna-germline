#!/usr/bin/env nextflow

include { CHECK_PRE_VALIDATION } from '../modules/utilities'
include { CONVERT_TO_CRAM } from '../modules/samtools'
include { VALIDATE_SAM_FILE } from '../modules/picard'

workflow BAM_TO_CRAM {
    take:
    sample_name
    bam_ch
    ref_fasta
    ref_fasta_index
    ref_dict
    duplication_metrics
    chimerism_metrics

    main:
    max_duplication_in_reasonable_sample = 0.30
    max_chimerism_in_reasonable_sample = 0.15

    cram_ch = CONVERT_TO_CRAM(bam_ch, ref_fasta, ref_fasta_index, sample_name)

    CHECK_PRE_VALIDATION(
        duplication_metrics,
        chimerism_metrics,
        max_duplication_in_reasonable_sample,
        max_chimerism_in_reasonable_sample,
    )

    is_outlier_data = CHECK_PRE_VALIDATION.out.is_outlier_data.text.trim().toBoolean()

    VALIDATE_SAM_FILE(
        cram_ch.cram,
        cram_ch.cram_index,
        ref_fasta,
        ref_fasta_index,
        ref_dict,
        1000000000,
        ["MISSING_TAG_NM"],
        is_outlier_data,
        "${sample_name}.cram.validation_report.txt",
    )

    emit:
    cram = cram_ch.cram
    cram_index = cram_ch.cram_index
    cram_md5 = cram_ch.cram_md5
    validation_report = VALIDATE_SAM_FILE.out.report
}
