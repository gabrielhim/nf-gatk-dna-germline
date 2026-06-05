#!/usr/bin/env nextflow

process APPLY_BQSR {

    container 'us.gcr.io/broad-gatk/gatk:4.3.0.0'
    memory 3500.MB

    input:
    path bam
    path bam_index
    path recalibration_report
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    val sequence_group_intervals
    val compression_level
    val bin_base_qualities
    val somatic
    val output_prefix

    output:
    path "${output_prefix}.recalibrated.${task.index}.bam", emit: recalibrated_bam
    path "${output_prefix}.recalibrated.${task.index}.bai", emit: recalibrated_bam_index
    path "${output_prefix}.recalibrated.${task.index}.bam.md5", emit: recalibrated_bam_md5

    script:
    def bin_base_qualities_arg = bin_base_qualities ? "--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30" : ""
    def bin_somatic_base_qualities_arg = bin_base_qualities && somatic ? "--static-quantized-quals 40 --static-quantized-quals 50" : ""
    def intervals_arg = sequence_group_intervals.collect { target -> "--intervals ${target}" }.join(' ')
    def java_initial_memory_mb = task.memory.toMega() - 500
    def java_max_memory_mb = task.memory.toMega() - 500
    """
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
        -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -Dsamjdk.compression_level=${compression_level} -Xms${java_initial_memory_mb}m -Xmx${java_max_memory_mb}m" \
        ApplyBQSR \
        --create-output-bam-md5 \
        --add-output-sam-program-record \
        --reference ${ref_fasta} \
        --input ${bam} \
        --use-original-qualities \
        --output ${output_prefix}.recalibrated.${task.index}.bam \
        --bqsr-recal-file ${recalibration_report} \
        ${bin_base_qualities_arg} \
        ${bin_somatic_base_qualities_arg} \
        ${intervals_arg}
    """
}

process BASE_RECALIBRATOR {

    container 'us.gcr.io/broad-gatk/gatk:4.3.0.0'
    memory 6000.MB

    input:
    path bam
    path bam_index
    path dbsnp_vcf
    path dbsnp_vcf_index
    path known_indels_sites_vcfs
    path known_indels_sites_indices
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    val sequence_group_intervals
    val output_prefix

    output:
    path "${output_prefix}.recal_data.${task.index}.txt", emit: recalibration_report

    script:
    def known_sites_arg = "--known-sites ${dbsnp_vcf} " + known_indels_sites_vcfs.collect { vcf -> "--known-sites ${vcf}" }.join(' ')
    def intervals_arg = sequence_group_intervals.collect { target -> "--intervals ${target}" }.join(' ')
    def java_initial_memory_mb = task.memory.toMega() - 1000
    def java_max_memory_mb = task.memory.toMega() - 500
    """
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
        -Xloggc:gc_log.log -Xms${java_initial_memory_mb}m -Xmx${java_max_memory_mb}m" \
        BaseRecalibrator \
        --reference ${ref_fasta} \
        --input ${bam} \
        --use-original-qualities \
        --output ${output_prefix}.recal_data.${task.index}.txt \
        ${known_sites_arg} \
        ${intervals_arg}
    """
}

process CALIBRATE_DRAGSTR_MODEL {

    container 'us.gcr.io/broad-gatk/gatk:4.3.0.0'
    memory 2500.MB
    cpus 4

    input:
    path alignment
    path alignment_index
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path str_table
    val output_prefix

    output:
    path "${output_prefix}.dragstr_model.txt", emit: dragstr_model

    script:
    def java_max_memory_mb = task.memory.toMega() < 2000 ? 1000 : task.memory.toMega() - 1000
    """
    gatk --java-options "-Xmx${java_max_memory_mb}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Dsamjdk.reference_fasta=${ref_fasta}" \
        CalibrateDragstrModel \
        --reference ${ref_fasta} \
        --input ${alignment} \
        --str-table-path ${str_table} \
        --output ${output_prefix}.dragstr_model.txt \
        --threads ${task.cpus}
    """
}

process DRAGEN_HARD_FILTER_VCF {

    container 'us.gcr.io/broad-gatk/gatk:4.6.1.0'
    memory 3000.MB

    input:
    path vcf
    path vcf_index
    val is_gvcf
    val output_prefix

    output:
    path "*.vcf.gz", emit: filtered_vcf
    path ".vcf.gz.tbi", emit: filtered_vcf_index

    script:
    def vcf_suffix = is_gvcf ? "hard-filtered.g.vcf.gz" : "hard-filtered.vcf.gz"
    """
    gatk --java-options "-Xms2000m -Xmx2500m" \
        VariantFiltration \
        --variant ${vcf} \
        --filter-expression "QUAL < 10.4139" \
        --filter-name "DRAGENHardQUAL" \
        --output ${output_prefix}.${vcf_suffix}
    """
}

process GATHER_BQSR_REPORTS {

    container 'us.gcr.io/broad-gatk/gatk:4.3.0.0'
    memory 3500.MB

    input:
    path bqsr_reports
    val output_prefix

    output:
    path "${output_prefix}.recal_data.txt", emit: bqsr_report

    script:
    def inputs_arg = bqsr_reports.collect {report -> " --input ${report}"}.join(' ')
    """
    gatk --java-options "-Xms3000m -Xmx3000m" \
        GatherBQSRReports \
        ${inputs_arg} \
        --output ${output_prefix}.recal_data.txt
    """
}

process HAPLOTYPE_CALLER {

    container 'us.gcr.io/broad-gatk/gatk:4.6.1.0'
    memory 8000.MB

    input:
    path(bam)
    path(bam_index)
    path(interval_list)
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path dragstr_model
    val make_gvcf
    val make_bamout
    val run_dragen_mode_variant_calling
    val use_spanning_event_genotyping
    val contamination
    val output_prefix

    output:
    path "*.vcf.gz", emit: vcf
    path "*.vcf.gz.tbi", emit: vcf_index
    path "*.bamout.bam", emit: bamout, optional: true

    script:
    def output_basename = "${output_prefix}.haplotype_caller.${task.index}"
    def vcf_suffix = make_gvcf ? "g.vcf.gz" : "vcf.gz"
    def dragen_mode_arg = run_dragen_mode_variant_calling ? "--dragen-mode" : ""
    def spanning_event_genotyping_arg = use_spanning_event_genotyping ? "" : "--disable-spanning-event-genotyping"
    def dragstr_model_arg = dragstr_model ? "--dragstr-params-path ${dragstr_model}" : ""
    def gvcf_gq_bands_arg = [10, 20, 30, 40, 50, 60, 70, 80, 90].collect {n -> "--gvcf-gq-bands ${n}"}.join(' ')
    def gvcf_arg = make_gvcf ? "--emit-ref-confidence GVCF --annotation-group AS_StandardAnnotation" : ""
    def bamout_arg = make_bamout ? "--bam-output ${output_prefix}.bamout.bam" : ""
    def java_memory_size_mb = task.memory.toMega() - 1000
    """
    gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        --reference ${ref_fasta} \
        --input ${bam} \
        --intervals ${interval_list} \
        --output ${output_basename}.${vcf_suffix} \
        --contamination-fraction-to-filter ${contamination} \
        --annotation-group StandardAnnotation \
        --annotation-group StandardHCAnnotation \
        ${dragen_mode_arg} \
        ${spanning_event_genotyping_arg} \
        ${dragstr_model_arg} \
        ${gvcf_gq_bands_arg} \
        ${gvcf_arg} \
        ${bamout_arg}
    """
}

process REBLOCK_GVCF {

    container 'us.gcr.io/broad-gatk/gatk:4.6.1.0'
    memory 3750.MB

    input:
    path gvcf
    path gvcf_index
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    val tree_score_cutoff
    val move_filters_to_genotypes
    val output_prefix

    output:
    path "${output_prefix}.reblocked.g.vcf.gz", emit: reblocked_gvcf
    path "${output_prefix}.reblocked.g.vcf.gz.tbi", emit: reblocked_gvcf_index

    script:
    def tree_score_cutoff_arg = tree_score_cutoff ? "--tree-score-threshold-to-no-call ${tree_score_cutoff}" : ""
    def move_filters_to_genotypes_arg = move_filters_to_genotypes ? "--add-site-filters-to-genotype" : ""
    """
    gatk --java-options "-Xms3000m -Xmx3000m" \
        ReblockGVCF \
        --reference ${ref_fasta} \
        --variant ${gvcf} \
        --do-qual-score-approximation \
        --floor-blocks \
        --gvcf-gq-bands 20 \
        --gvcf-gq-bands 30 \
        --gvcf-gq-bands 40 \
        ${tree_score_cutoff_arg} \
        ${move_filters_to_genotypes_arg} \
        --output ${output_prefix}.reblocked.g.vcf.gz
    """
}

process VALIDATE_VCF {

    container 'us.gcr.io/broad-gatk/gatk:4.6.1.0'
    memory 7000.MB

    input:
    path vcf
    path vcf_index
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path dbsnp_vcf
    path dbsnp_vcf_index
    path calling_interval_list
    val is_gvcf
    val extra_args

    script:
    def dbsnp_arg = dbsnp_vcf ? "--dbsnp ${dbsnp_vcf}" : ""
    def gvcf_arg = is_gvcf ? "--validate-GVCF" : ""
    def java_memory_mb = task.memory.toMega() - 2000
    """
    gatk --java-options "-Xms${java_memory_mb}m -Xmx${java_memory_mb}m" \
        ValidateVariants \
        --variant ${vcf} \
        --reference ${ref_fasta} \
        --intervals ${calling_interval_list} \
        ${gvcf_arg} \
        --validation-type-to-exclude ALLELES \
        ${dbsnp_arg} \
        ${extra_args}
    """
}
