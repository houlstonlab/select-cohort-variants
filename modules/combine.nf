process COMBINE {
    tag "${pheno}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/combined/", mode: 'copy')

    input:
    tuple val(pheno), val(chrom), val(category),
          path(file), path(index), env(n_vars)

    output:
    tuple val(pheno), val(category),
          path("${pheno}.${category}.vcf.gz"),
          path("${pheno}.${category}.vcf.gz.tbi"),
		  path("${pheno}.${category}.annotations.tsv"),
          env(n_vars)
        
    script:
    """
    #!/bin/bash
    # Combine vcfs and update IDs
    bcftools concat --naive ${file} | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
    bcftools view --threads ${task.cpu} -Oz -o ${pheno}.${category}.vcf.gz
    
    # Index the vcf
    tabix ${pheno}.${category}.vcf.gz

	# Count the number of variants
    n_vars=\$(bcftools index -n ${pheno}.${category}.vcf.gz)

	# Extract annotations
	bcftools +split-vep \
		-s worst \
		-c Gene \
		-f '%CHROM:%POS:%REF:%ALT\t%Gene\t%CSQ\n' \
		-d -A tab \
		${pheno}.${category}.vcf.gz \
		> ${pheno}.${category}.annotations.tsv
	"""
}
