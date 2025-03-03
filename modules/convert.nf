process CONVERT {
    tag "${pheno}:${category}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/plinked", mode: 'copy')

    input:
    tuple val(pheno), val(category),
          path(file), path(index), 
          path(annotations), env(n_vars)

    output:
    tuple val(pheno), val(category),
          path("${pheno}.${category}.bim"),
          path("${pheno}.${category}.bed"),
          path("${pheno}.${category}.fam"),
          path("${pheno}.${category}.nosex"),
          path("${pheno}.${category}.log"),
          path(annotations)

    script:
    """
    #!/bin/bash
    plink \
        --vcf ${file} \
        --make-bed \
        --const-fid 0 \
        --out ${pheno}.${category}
    """
}