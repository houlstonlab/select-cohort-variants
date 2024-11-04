process COORDINATES {
    tag "${chrom}"

    label 'simple'

    container  params.bioconductor

    publishDir("${params.output_dir}/coordiantes", mode: 'copy')

    input:
    tuple val(chrom), val(genome), val(style)

    output:
    tuple val(chrom), path("${chrom}.bed")
 
    script:
    """
    #!/bin/bash
    generate_coordinates.R ${genome} ${style} ${chrom} ${chrom}.bed
    """
}
