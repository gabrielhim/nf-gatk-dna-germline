#!/usr/bin/env nextflow

process BWA_MEM_AND_MERGE_ALIGNMENT {

    container 'us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748'
    memory 14.GB
    cpus 16

    input:
    tuple path(unmapped_bam), val(output_prefix)
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path ref_alt
    path ref_amb
    path ref_ann
    path ref_bwt
    path ref_pac
    path ref_sa

    val compression_level
    val allow_empty_ref_alt
    val hard_clip_reads
    val unmap_contaminant_reads

    output:
    path "${output_prefix}.merged.bam", emit: merged_bam
    path "${output_prefix}.bwa.stderr.log", emit: bwa_stderr_log

    script:
    def bwa_command_line = "bwa mem -K 100000000 -p -v 3 -t 16 -Y \$bash_ref_fasta"
    def hard_clip_reads_arg = hard_clip_reads ? "CLIP_OVERLAPPING_READS=true CLIP_OVERLAPPING_READS_OPERATOR=H'" : ""
    """
    # This is done before "set -o pipefail" because "bwa" will have a rc=1 and we don't want to allow rc=1 to succeed
    # because the sed may also fail with that error and that is something we actually want to fail on.
    BWA_VERSION=\$(/usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //')

    set -o pipefail
    set -e

    if [ -z \${BWA_VERSION} ]; then
        exit 1;
    fi

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it or allow_empty_ref_alt is set
    if [ -s ${ref_alt} ] || ${allow_empty_ref_alt}; then
        java -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
            SamToFastq \
            INPUT=${unmapped_bam} \
            FASTQ=/dev/stdout \
            INTERLEAVE=true \
            NON_PF=true | \
        /usr/gitc/${bwa_command_line} /dev/stdin - 2> >(tee ${output_prefix}.bwa.stderr.log >&2) | \
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
            PROGRAM_RECORD_ID="bwamem" \
            PROGRAM_GROUP_VERSION="\${BWA_VERSION}" \
            PROGRAM_GROUP_COMMAND_LINE="${bwa_command_line}" \
            PROGRAM_GROUP_NAME="bwamem" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=${unmap_contaminant_reads} \
            ADD_PG_TAG_TO_READS=false

      if ${!allow_empty_ref_alt}; then
            grep -m1 "read .* ALT contigs" ${output_prefix}.bwa.stderr.log | \
            grep -v "read 0 ALT contigs"
      fi

    # else ref_alt is empty or could not be found
    else
        echo ref_alt input is empty or not provided. >&2
        exit 1;
    fi
    """
}
