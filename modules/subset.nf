process SUBSET {
    tag "${pheno}:${chrom}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/pheno", mode: 'symlink')

    input:
    tuple val(pheno), val(chrom), path(gene_coords),
          path(file), path(index),
          path(cases)

    output:
    tuple val(pheno), val(chrom),
          path("${pheno}.${chrom}.vcf.gz"),
          path("${pheno}.${chrom}.vcf.gz.tbi")

    script:
    """
    #!/bin/bash
    # Subset pheno
    bcftools view \
        -S ${cases} \
        -R ${gene_coords} \
        -i 'FILTER="PASS"' \
        -g het \
        ${file} | \
    bcftools norm -m -any | \
    bcftools +fill-tags -- -t all | \
    bcftools +setGT -- -t . -n 0 | \
    bcftools view -g het --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.vcf.gz

    tabix ${pheno}.${chrom}.vcf.gz
    """
}
