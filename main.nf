#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILL }        from './modules/fill.nf'
include { FILTER }      from './modules/filter.nf'
include { COMBINE }     from './modules/combine.nf'
include { CONVERT }     from './modules/convert.nf'
include { EXTRACT }     from './modules/extract.nf'
include { AGGREGATE }   from './modules/aggregate.nf'

// Define input channels
variants_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [
        row.cohort,
        row.chrom ?: (1..22).collect { "chr$it" } + ['chrX', 'chrY'],
        file(row.file), file(row.index),
        file(row.samples)
    ] }
    | transpose

genes_coords_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ 
        row.cohort,
        row.chrom ?: (1..22).collect { "chr$it" } + ['chrX', 'chrY'],
        params.genome, params.style
    ] }
    | transpose
    | groupTuple(by: [1,2,3])

category_ch = Channel.of(params.categories.split(','))
variable_ch = Channel.of( 'rlist', 'snplist', 'frqx' )
// variable_ch = Channel.of( 'annotations', 'list', 'rlist', 'snplist', 'frqx' )

workflow  {
    // Subset, Filter and Extract qualifying variants
    genes_coords_ch
        | COORDINATES
        | transpose
        | combine(variants_ch, by: [0,1])
        | SUBSET
        | ( params.fill ? FILL : map {it} )
        | combine(category_ch)
        | FILTER
        | filter { it[5].toInteger() > 0 }
        | groupTuple(by: [0,2])
        | COMBINE
        | CONVERT
        | combine(variable_ch)
        | EXTRACT
        | filter { it[3] == 'rlist' } 
        | AGGREGATE
}
