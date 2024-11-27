process COMBINE {
    tag "${pheno}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/combined/", mode: 'copy')

    input:
    tuple val(pheno), val(chrom), val(category),
          path(file), path(index), path(variants)

    output:
    tuple val(pheno), val(category),
          path("${pheno}.${category}.vcf.gz"),
          path("${pheno}.${category}.vcf.gz.tbi"),
          path("${pheno}.${category}.tsv")
        
    script:
    """
    #!/bin/bash
    # Combine vcfs
    bcftools concat \
        --naive \
        ${file} \
        --threads ${task.cpu} \
        -Oz -o ${pheno}.${category}.vcf.gz
    tabix ${pheno}.${category}.vcf.gz

    # Combine variants
    cat ${variants} > ${pheno}.${category}.tsv
    """
}
