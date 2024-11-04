process REPORT {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/reports", mode: 'copy')

    input:
    tuple val(pheno), val(chrom), val(category),
          path(file), path(index), path(variants)

    output:
    tuple val(pheno), val(chrom), val(category), val('report'),
          path("${pheno}.${chrom}.${category}.report.txt")

    script:
    """
    #!/bin/bash
    calculate_stats.sh \
        ${pheno} ${chrom} ${category} \
        ${file} \
        > ${pheno}.${chrom}.${category}.report.txt
    """
}
