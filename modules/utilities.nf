#!/usr/bin/env nextflow

process CHECK_PRE_VALIDATION {

    container 'python:3.12-slim'
    memory 2.GB

    input:
    path duplication_metrics
    path chimerism_metrics
    val max_duplication_in_reasonable_sample
    val max_chimerism_in_reasonable_sample

    output:
    path "duplication_value.txt", emit: duplication_rate
    path "chimerism_value.txt", emit: chimerism_rate
    path "is_outlier_data.txt", emit: is_outlier_data

    script:
    """
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ${duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ${chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python3 <<CODE
    import csv
    with open("duplication.csv") as dupfile:
        reader = csv.DictReader(dupfile, delimiter="\\t")
        for row in reader:
            with open("duplication_value.txt", "w") as file:
                file.write(row["PERCENT_DUPLICATION"])

    with open("chimerism.csv") as chimfile:
        reader = csv.DictReader(chimfile, delimiter="\\t")
        for row in reader:
            with open("chimerism_value.txt", "w") as file:
                file.write(row["PCT_CHIMERAS"])
    CODE

    if (( \$(cat "duplication_value.txt") > ${max_duplication_in_reasonable_sample} )) || \
        (( \$(cat "chimerism_value.txt") > ${max_chimerism_in_reasonable_sample} )); then
        echo "true" > is_outlier_data.txt
    else
        echo "false" > is_outlier_data.txt
    fi
    """
}

process CREATE_SEQUENCE_GROUPING_TSV {

    container 'python:3.12-slim'
    memory 2.GB

    input:
    path ref_dict

    output:
    path 'sequence_grouping.txt', emit: sequence_grouping
    path 'sequence_grouping_with_unmapped.txt', emit: sequence_grouping_with_unmapped

    script:
    """
    #!/usr/bin/env python3

    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.

    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt", "w") as tsv_file:
        tsv_file.write(tsv_string)

    tsv_string += "\\n" + "unmapped"
    with open("sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
        tsv_file_with_unmapped.write(tsv_string)
    """
}
