#!/usr/bin/env nextflow

process CHECK_CONTAMINATION {

    container 'us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840'
    memory '7.5 GB'

    input:
    path bam
    path bam_index
    path contamination_sites_ud
    path contamination_sites_bed
    path contamination_sites_mu
    path ref_fasta
    path ref_fasta_index
    val contamination_underestimation_factor
    val compression_level
    val disable_sanity_check
    val output_prefix

    output:
    path "${output_prefix}.selfSM", emit: per_sample_stats
    stdout emit: contamination

    script:
    def disable_sanity_check_arg = disable_sanity_check ? "--DisableSanityCheck" : ""
    """
    set -e

    # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
        --Verbose \
        --NumPC 4 \
        --Output ${output_prefix} \
        --BamFile ${bam} \
        --Reference ${ref_fasta} \
        --UDPath ${contamination_sites_ud} \
        --MeanPath ${contamination_sites_mu} \
        --BedPath ${contamination_sites_bed} \
        ${disable_sanity_check_arg} \
        1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open("${output_prefix}.selfSM") as selfSM:
        reader = csv.DictReader(selfSM, delimiter="\\t")
        i = 0
        for row in reader:
            if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
                # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
                # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
                # vcf and bam.
                sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
                sys.exit(1)
        print(float(row["FREEMIX"])/${contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
            sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
            sys.exit(2)
    CODE
    """
}