process SUBSET {
    tag "${pheno}:${chrom}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/cohort", mode: 'symlink')

    input:
    tuple val(chrom), path(file), path(index),
          path(gene_coords)

    output:
    tuple val(chrom),
          path("gnomad.${chrom}.vcf.gz"),
          path("gnomad.${chrom}.vcf.gz.tbi")

    script:
    """
    #!/bin/bash
    # Subset cohort
    bcftools view \
        -S ${pheno_file} \
        -R ${gene_coords} \
        -i 'FILTER="PASS"' \
        -g het \
        ${cohort_file} | \
    bcftools norm -m -any | \
    bcftools +fill-tags -- -t all | \
    bcftools +setGT -- -t . -n 0 | \
    bcftools +setGT -- -t q -n 0 -i 'FMT/GQ < ${params.GQ} | FMT/DP < ${params.DP} | VAF < ${params.VAF}' | \
    bcftools +fill-tags -- -t all | \
    bcftools view -g het --threads ${task.cpu} -Oz -o gnomad.${chrom}.vcf.gz

    tabix gnomad.${chrom}.vcf.gz
    """
}
