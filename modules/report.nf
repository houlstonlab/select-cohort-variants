process REPORT {
    tag "${pheno}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/reports", mode: 'copy')

    input:
    tuple val(pheno), val(category),
          path(file), path(index), path(variants)

    output:
    tuple val(pheno), val(category), val('report'),
          path("${pheno}.${category}.report.txt")

    script:
    """
    #!/bin/bash
    calculate_stats.sh \
        ${pheno} ${category} ${file} \
        > ${pheno}.${category}.report.txt
    """
}
