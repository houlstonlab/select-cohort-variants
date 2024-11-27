process EXTRACT {
    tag "${pheno}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/variants", mode: 'copy')

    input:
    tuple val(pheno), val(category),
          path(file), path(index), path(variants),
          val(variable)

    output:
    tuple val(pheno), val(category), val(variable),
	      path("${pheno}.${category}.${variable}.tsv")

    script:
    if ( variable == 'variants' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\n' \
			${file} \
			> "${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'annotations' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%CSQ\n' \
			-d -A tab \
			${file} \
			> "${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'genotypes' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]\n' \
			${file} \
			> "${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'frequency' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]\n' \
			${file} | \
			summarize_genotypes.awk \
			> "${pheno}.${category}.${variable}.tsv"
		"""
    } 
}
