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
    path "${output_prefix}.recal_data.${task.index}.csv", emit: recalibration_report

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
        --output ${output_prefix}.recal_data.${task.index}.csv \
        ${known_sites_arg} \
        ${intervals_arg}
    """
}

process GATHER_BQSR_REPORTS {

    container 'us.gcr.io/broad-gatk/gatk:4.3.0.0'
    memory 3500.MB

    input:
    path bqsr_reports
    val output_prefix

    output:
    path "${output_prefix}.recal_data.csv", emit: bqsr_report

    script:
    def inputs_arg = bqsr_reports.collect {report -> " --input ${report}"}.join(' ')
    """
    gatk --java-options "-Xms3000m -Xmx3000m" \
        GatherBQSRReports \
        ${inputs_arg} \
        --output ${output_prefix}.recal_data.csv
    """
}
