#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(readr)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 2) {
    stop("Usage: Rscript 2_table_to_matrix.R <input_file.txt> <output_file.txt>")
  }

  input_file <- argv[1]
  output_file <- argv[2]

  # Load data
  dat <- read.delim(input_file, stringsAsFactors = FALSE)

  # Pivot into wide format
  sorted <- dat[order(dat$Sample_Id), ]
  mat <- pivot_wider(sorted, names_from = Sample_Id, values_from = Count)

  # Write to .txt file
  write.table(mat, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

main()

