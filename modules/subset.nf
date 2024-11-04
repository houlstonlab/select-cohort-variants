process SUBSET {
    tag "${pheno}:${chrom}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/cohort", mode: 'symlink')

    input:
    tuple val(chrom), path(gene_coords),
          val(cohort), path(cohort_file), path(cohort_index),
          val(pheno), path(pheno_file)

    output:
    tuple val(pheno), val(chrom),
          path("${pheno}.${chrom}.vcf.gz"),
          path("${pheno}.${chrom}.vcf.gz.tbi")

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
    bcftools view -g het --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.vcf.gz

    tabix ${pheno}.${chrom}.vcf.gz
    """
}
