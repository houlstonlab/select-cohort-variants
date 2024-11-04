#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
style <- args[2]
chrom <- args[3]
chrom.coordiantes <- args[4]

# load genes
# genome <- 'hg38'
# chrom <- 'chr21'
# chrom.coordiantes <- paste0('tmp/', chrom, '.bed')

if (genome == 'hg38') txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# elseif (genome == 'hg37') txdb <- TxDb.Hsapiens.UCSC.hg37.knownGene::TxDb.Hsapiens.UCSC.hg37.knownGene

gene_coordinates <- GenomicFeatures::genes(
  txdb,
  filter = list(tx_chrom = chrom),
  columns = AnnotationDbi::columns(txdb)
)

GenomeInfoDb::seqlevels(gene_coordinates) <- chrom
GenomeInfoDb::seqlevelsStyle(gene_coordinates) <- style

rtracklayer::export.bed(
  gene_coordinates,
  chrom.coordiantes
)
