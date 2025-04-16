#!/usr/bin/env Rscript

# Load only required libraries
library(sva)
library(readr)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 2) {
    stop("Usage: Rscript combat_seq_correct.R <input_file.txt> <output_file.txt>")
  }

  input_path <- argv[1]
  output_path <- argv[2]

  # Read input count table
  df <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

  # Step 1: Rename specific 0_CAD samples to 00_CAD (for batch assignment)
  samples_to_rename <- c("0_CAD_KGKONTR16", "0_CAD_KGKONTR17", "0_CAD_KGKONTR23",
                         "0_CAD_KGKONTR39", "0_CAD_KGKONTR42")
  colnames(df) <- ifelse(colnames(df) %in% samples_to_rename,
                         sub("^0_CAD", "00_CAD", colnames(df)),
                         colnames(df))

  # Step 2: Set miRNA names as row names and drop miRNA column
  rownames(df) <- df$miRNA
  df <- df[, -1]

  # Step 3: Assign batch labels
  batch <- ifelse(grepl("^00_CAD|^1_CAD|^2_CAD", colnames(df)), 1, 2)

  # Step 4: Apply ComBat-seq
  adjusted_counts <- ComBat_seq(as.matrix(df), batch = batch, group = NULL)

  # Step 5: Rename 00_CAD back to 0_CAD
  colnames(adjusted_counts) <- sub("^00_CAD", "0_CAD", colnames(adjusted_counts))

  # Step 6: Write output with miRNA column
  adjusted_df <- data.frame(miRNA = rownames(adjusted_counts), adjusted_counts, check.names = FALSE)
  write.table(adjusted_df, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

  message("✅ ComBat-seq correction complete. Output saved to: ", output_path)
}

main()

