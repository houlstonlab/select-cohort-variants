#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { EXTRACT }     from './modules/extract.nf'
include { AGGREGATE }   from './modules/aggregate.nf'
include { REPORT }      from './modules/report.nf'

// Define input channels
genes_coords_ch =  Channel.of (1..22, 'X', 'Y')
    | map { [ "chr${it}", params.genome, params.style ] }

variants_ch = Channel.fromFilePairs(params.vcf, flat: true)

cases_ch = Channel.fromPath(params.cases)
    | map { [it.simpleName, it] }

category_ch = Channel.of( 'Pathogenic', 'Damaging', 'Splicing',
                        //   'NotScored', 'Scored',
                          'High', 'PTV', 'Stop' )
variable_ch = Channel.of( 'variants', 'annotations', 'genotypes', 'frequency' )

workflow  {
    // Combine genes, variants, and cases
    // Subset, Filter and Extract qualifying variants
    genes_coords_ch
        | COORDINATES

    variants_ch
        | combine(cases_ch, by: 0)
        | combine(COORDINATES.out)
        | SUBSET
        | combine(category_ch)
        | FILTER
        | combine(variable_ch)
        | EXTRACT
        | filter { it[3] == 'frequency' }
        | multiMap {
            cat: it
            all: [it[0], it[1], 'ALL', it[3], it[4]]
        } 
        | set { variants }

    // Concatenate and group genotypes by category or 'all'
    variants.cat 
        | concat(variants.all) 
        | groupTuple(by: [0,1,2,3])
        | AGGREGATE

    // Stats
    FILTER.out | REPORT

    // Collect and store the summary files
    EXTRACT.out
        | concat(AGGREGATE.out)
        | concat(REPORT.out)
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/summary",
        )
        { item -> [ "${item[0]}.${item[2]}.${item[3]}.tsv", item.last() ] } 
}
