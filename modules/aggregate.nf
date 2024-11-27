process AGGREGATE {
    tag "${pheno}:${category}"

    label 'simple'

    container params.rocker

    publishDir("${params.output_dir}/aggregate", mode: 'copy')

    input:
    tuple val(pheno), val(category), val(variable),
          path(file)

    output:
    tuple val(pheno), val(category), val('aggregate'),
          path("${pheno}.${category}.aggregate.tsv")
 
    script:
    """
    #!/bin/bash
	cat ${file} | awk '!seen[\$0]++' > variants.txt
    aggregate_genotyeps.R variants.txt ${pheno}.${category}.aggregate.tsv
    """
}