#!/usr/bin/env Rscript

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Function to count cases
count_cases <- function(x, compound = FALSE) {
  # collapse all cases to a string
  s <- paste(na.omit(x), collapse = ',')
  
  if (all(is.na(s))) {
    # if no cases, return 0
    count <- 0
  } else {
    if ( compound ) {
      # cases showing up more than once, is compound
      count <- sum(table(unlist(strsplit(s, ','))) > 1)
    } else {
      # the number of unique cases
      count <- length(unique(unlist(strsplit(s, ','))))
    }
  }
  return(count)
}

# Read input file
if (file.size(input_file) > 0) {
  d <- readr::read_tsv(input_file, col_names = FALSE)
  
  # Process data
  d2 <- dplyr::group_by(d, gene = X1)
  d2 <- dplyr::summarise(d2, het = count_cases(X4),
                  hom = count_cases(X6),
                  ch  = count_cases(X4, compound = TRUE))
} else {
  d2 <- dplyr::tibble(gene = character(), het = numeric(), hom = numeric(), ch = numeric())
}

# Write output file
readr::write_tsv(d2, output_file)
