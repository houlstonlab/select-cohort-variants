process AGGREGATE {
    tag "${chrom}:${category}"

    label 'simple'

    container params.rocker

    publishDir("${params.output_dir}/aggregate", mode: 'copy')

    input:
    tuple val(pheno), val(chrom), val(category), val(variable),
          path(file)

    output:
    tuple val(pheno), val(chrom), val(category), val('aggregate'),
          path("${pheno}.${chrom}.${category}.aggregate.tsv")
 
    script:
    """
    #!/bin/bash
	cat ${file} | awk '!seen[\$0]++' > variants.txt
    aggregate_genotyeps.R variants.txt ${pheno}.${chrom}.${category}.aggregate.tsv
    """
}