#!/usr/bin/env nextflow

process DRAGMAP_AND_MERGE_ALIGNMENT {

    container 'us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530'
    memory 40.GB
    cpus 16

    input:
    tuple val(output_prefix), path(unmapped_bam)
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path reference_bin
    path hash_table_cfg_bin
    path hash_table_cmp

    val compression_level
    val hard_clip_reads
    val unmap_contaminant_reads

    output:
    path "${output_prefix}.merged.bam", emit: merged_bam
    path "${output_prefix}.dragmap.stderr.log", emit: dragmap_stderr_log

    script:
    def dragmap_command_line = "dragen-os -b \$input_bam -r dragen_reference --interleaved=1 --preserve-map-align-order true"
    def hard_clip_reads_arg = hard_clip_reads ? "CLIP_OVERLAPPING_READS=true CLIP_OVERLAPPING_READS_OPERATOR=H'" : ""
    """
   DRAGMAP_VERSION=\$(dragen-os --version)

    if [ -z \${DRAGMAP_VERSION} ]; then
        exit 1;
    fi

    mkdir dragen_reference
    mv ${reference_bin} ${hash_table_cfg_bin} ${hash_table_cmp} dragen_reference

    input_bam=${unmapped_bam}
    ${dragmap_command_line} 2> >(tee ${output_prefix}.dragmap.stderr.log >&2) | samtools view -h -O BAM - > aligned.bam

    java -Dsamjdk.compression_level=${compression_level} -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${unmapped_bam} \
        OUTPUT=${output_prefix}.merged.bam \
        REFERENCE_SEQUENCE=${ref_fasta} \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        ${hard_clip_reads_arg} \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="dragen-os" \
        PROGRAM_GROUP_VERSION="\${DRAGMAP_VERSION}" \
        PROGRAM_GROUP_COMMAND_LINE="${dragmap_command_line}" \
        PROGRAM_GROUP_NAME="dragen-os" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=${unmap_contaminant_reads} \
        ADD_PG_TAG_TO_READS=false
    """
}
