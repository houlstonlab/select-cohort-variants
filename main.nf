#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { COMBINE }     from './modules/combine.nf'
include { EXTRACT }     from './modules/extract.nf'
include { AGGREGATE }   from './modules/aggregate.nf'
include { REPORT }      from './modules/report.nf'

// Define input channels
variants_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ row.cohort, file(row.file), file(row.index), file(row.samples) ] }

genes_coords_ch = Channel.fromPath(params.cohorts)
    | splitCsv(header: true, sep: ',')
    | map { row -> [ 
        row.chrom ?: (1..22).collect { "chr$it" } + ['chrX', 'chrY'],
        params.genome, params.style
    ] }
    | transpose()

category_ch = Channel.of(params.categories.split(','))
variable_ch = Channel.of( 'variants', 'annotations', 'genotypes', 'frequency' )

workflow  {
    // Combine genes, variants, and cases
    // Subset, Filter and Extract qualifying variants
    genes_coords_ch
        | COORDINATES

    variants_ch
        | combine(COORDINATES.out)
        | SUBSET
        | combine(category_ch)
        | FILTER
        | groupTuple(by: [0,2])
        | COMBINE
        | combine(variable_ch)
        | EXTRACT
        | filter { it[2] == 'frequency' }
        | multiMap {
            cat: it
            all: [it[0], 'ALL', it[2], it[3]]
        } 
        | set { variants }

    // Concatenate and group genotypes by category or 'all'
    variants.cat 
        | concat(variants.all) 
        | groupTuple(by: [0,1,2])
        | AGGREGATE

    // Stats
    COMBINE.out | REPORT
}
