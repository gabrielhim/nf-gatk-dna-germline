#!/usr/bin/env nextflow

process SUBSET_CONTAMINATION_RESOURCES {

    container 'us.gcr.io/broad-gotc-prod/bedtools:2.27.1'
    memory '3.5 GB'

    input:
    path target_interval_list
    path contamination_sites_ud
    path contamination_sites_bed
    path contamination_sites_mu

    output:
    path "${contamination_sites_ud.baseName}.target.UD", emit: subset_ud
    path "${contamination_sites_bed.baseName}.target.bed", emit: subset_bed
    path "${contamination_sites_mu.baseName}.target.mu", emit: subset_mu
    path "target_overlap_counts.txt", emit: target_overlap_counts

    script:
    def output_ud = "${contamination_sites_ud.baseName}.target.UD"
    def output_bed = "${contamination_sites_bed.baseName}.target.bed"
    def output_mu = "${contamination_sites_mu.baseName}.target.mu"
    def target_overlap_counts = "target_overlap_counts.txt"
    """
    set -e -o pipefail

    grep -vE "^@" ${target_interval_list} |
       awk -v OFS='\t' '\$2=\$2-1' |
       /app/bedtools intersect -c -a ${contamination_sites_bed} -b - |
       cut -f6 > ${target_overlap_counts}

    function restrict_to_overlaps() {
        # print lines from whole-genome file from loci with non-zero overlap
        # with target intervals
        wgs_file=\$1
        exome_file=\$2
        paste ${target_overlap_counts} \$wgs_file |
            grep -Ev "^0" |
            cut -f 2- > \$exome_file
        echo "Generated \$exome_file"
    }

    restrict_to_overlaps ${contamination_sites_ud} ${output_ud}
    restrict_to_overlaps ${contamination_sites_bed} ${output_bed}
    restrict_to_overlaps ${contamination_sites_mu} ${output_mu}
    """
}
