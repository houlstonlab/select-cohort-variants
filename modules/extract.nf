process EXTRACT {
    tag "${pheno}:${category}:${variable}"

    label 'simple'

    container params.plink

    publishDir("${params.output_dir}/variants", mode: 'copy')

    input:
    tuple val(pheno), val(category),
          path(bim), path(bed), path(fam), path(nosex), path(log),
		  path(annotations), 
		  val(variable)

    output:
    tuple val(pheno), val(category), path(annotations), val(variable),
	      path("${pheno}.${category}.extracted.${variable}"),
	      path("${pheno}.${category}.extracted.log")

    script:
    if ( variable == 'list' ) {
		"""
		#!/bin/bash
		plink --bfile ${bim.baseName} --recode list --out ${pheno}.${category}.extracted
		"""
    } else if ( variable == 'rlist' ) {
		"""
		#!/bin/bash
		plink --bfile ${bim.baseName} --recode rlist --out ${pheno}.${category}.extracted
		"""
    } else if ( variable == 'snplist' ) {
		"""
		#!/bin/bash
		plink --bfile ${bim.baseName} --write-snplist --out ${pheno}.${category}.extracted
		"""
    } else if ( variable == 'frqx' ) {
		"""
		#!/bin/bash
		plink --bfile ${bim.baseName} --freqx --out ${pheno}.${category}.extracted
		"""
    } else {
		println "Variable ${variable} not recognized"
	}
}
