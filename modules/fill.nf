process FILL {
    tag "${pheno}:${chrom}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/filled", mode: 'symlink')

    input:
    tuple val(pheno), val(chrom),
          path(file), path(index)

    output:
    tuple val(pheno), val(chrom),
          path("${pheno}.${chrom}.filled.vcf.gz"),
          path("${pheno}.${chrom}.filled.vcf.gz.tbi")

    script:
    """
    #!/bin/bash
    # Subset pheno
    bcftools view ${file} | \
    bcftools +setGT -- -t q -n 0 -i 'FMT/GQ < ${params.GQ} | FMT/DP < ${params.DP} | VAF < ${params.VAF}' | \
    bcftools +fill-tags -- -t all | \
    bcftools view -g het --threads ${task.cpu} -Oz -o ${pheno}.${chrom}.filled.vcf.gz

    tabix ${pheno}.${chrom}.filled.vcf.gz
    """
}
