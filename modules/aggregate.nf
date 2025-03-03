process AGGREGATE {
    tag "${pheno}:${category}"

    label 'simple'

    container params.rocker

    publishDir("${params.output_dir}/aggregate", mode: 'copy')

    input:
    tuple val(pheno), val(category), path(annotations), 
          val(rlist), path(rlist_file), path(rlist_log)

    output:
    tuple val(pheno), val(category), 
          path("${pheno}.${category}.aggregate.tsv")
 
    script:
    """
    #!/bin/bash
	aggregate_genotyeps.R ${annotations} ${rlist_file} ${pheno}.${category}.aggregate.tsv
    """
}