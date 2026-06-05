#!/usr/bin/env nextflow

process CONVERT_TO_CRAM {

    container 'us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616'
    memory 3.GB

    input:
    path bam
    path ref_fasta
    path ref_fasta_index
    val output_prefix

    output:
    path "${output_prefix}.cram", emit: cram
    path "${output_prefix}.cram.crai", emit: cram_index
    path "${output_prefix}.cram.md5", emit: cram_md5

    script:
    """
    set -e
    set -o pipefail

    samtools view -C -T ${ref_fasta} ${bam} | \
    tee ${output_prefix}.cram | \
    md5sum | awk '{print \$1}' > ${output_prefix}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ${ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ${output_prefix}.cram
    """
}

// Use samtools to merge bamouts because using Picard produces an error.
process MERGE_BAMOUTS {

    container 'biocontainers/samtools:1.3.1'
    memory 4.GB

    input:
    path bams
    val output_filename

    output:
    path output_filename, emit: merged_bamout
    path "${output_filename}.bai", emit: merged_bamout_index

    script:
    """
    set -eo pipefail

    # Order files numerically using the indices in the filename to get a sorted output.
    bash_inputs_arg=\$(ls ${bams.join(' ')} | sort -V | tr '\\n' ' ')

    samtools merge ${output_filename} \$bash_inputs_arg
    samtools index ${output_filename}
    """
}
