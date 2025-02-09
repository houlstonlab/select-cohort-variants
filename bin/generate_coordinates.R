#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
genome            <- args[1]
style             <- args[2]
chrom             <- args[3]
coding            <- args[4]
chrom.coordiantes <- args[5]

# load genes
if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# elseif (genome == 'hg37') txdb <- TxDb.Hsapiens.UCSC.hg37.knownGene::TxDb.Hsapiens.UCSC.hg37.knownGene

if ( coding == 'true' ) {
  gene_coordinates <- GenomicFeatures::genes(
    txdb,
    filter = list(tx_chrom = chrom),
    columns = AnnotationDbi::columns(txdb)
  )
} else if ( coding == 'false' ) {
  gene_coordinates <- GenomicRanges::GRanges(
    seqnames = chrom,
    IRanges::IRanges(1, GenomeInfoDb::seqlengths(txdb)[chrom])
  )
} else {
  stop("coding can be 'true' or 'false'.")
}

GenomeInfoDb::seqlevels(gene_coordinates) <- chrom
GenomeInfoDb::seqlevelsStyle(gene_coordinates) <- style

rtracklayer::export.bed(
  gene_coordinates,
  chrom.coordiantes
)
