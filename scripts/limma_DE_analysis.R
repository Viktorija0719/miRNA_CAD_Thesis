#!/usr/bin/env Rscript

library(limma)
library(edgeR)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 2) {
    stop("Usage: Rscript limma_DE_analysis.R <input_file.txt> <output_directory>")
  }
  
  input_file <- argv[1]
  output_dir <- argv[2]
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load normalized and batch-corrected logCPM counts
  countTable <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t", header = TRUE)
  
  # Extract sample groups by prefix
  col_names <- colnames(countTable)
  group1_cols <- grep("^0_CAD", col_names)
  group2_cols <- grep("^1_CAD", col_names)
  group3_cols <- grep("^2_CAD", col_names)
  
  counts <- countTable[, c(group1_cols, group2_cols, group3_cols)]
  rownames(counts) <- countTable$miRNA  # Ensure rownames match miRNAs
  
  # Redefine indices after combining
  group1_idx <- 1:length(group1_cols)
  group2_idx <- (max(group1_idx) + 1):(max(group1_idx) + length(group2_cols))
  group3_idx <- (max(group2_idx) + 1):(max(group2_idx) + length(group3_cols))
  
  # Group assignment
  groups <- as.factor(c(
    rep("Group0", length(group1_idx)),
    rep("Group1", length(group2_idx)),
    rep("Group2", length(group3_idx))
  ))
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  
  logCPM <- as.matrix(counts)
  
  # Median expression
  medianGroup0 <- apply(logCPM[, group1_idx], 1, median)
  medianGroup1 <- apply(logCPM[, group2_idx], 1, median)
  medianGroup2 <- apply(logCPM[, group3_idx], 1, median)
  medianAll <- apply(logCPM, 1, median)
  
  # Filter
  selected <- which(medianAll > 3.0)
  logCPMFiltered <- logCPM[selected, ]
  
  # Linear modeling
  fit <- lmFit(logCPMFiltered, design)
  fit$genes <- rownames(logCPMFiltered)
  
  contrast_matrix <- makeContrasts(
    Group1vsGroup2 = Group1 - Group2,
    Group2vsGroup0 = Group2 - Group0,
    Group1vsGroup0 = Group1 - Group0,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  eBayesFit <- eBayes(fit2, trend = TRUE)
  
  comparisons <- c("Group1vsGroup2", "Group2vsGroup0", "Group1vsGroup0")
for (comp in comparisons) {
  res <- topTable(eBayesFit, coef = comp, number = Inf, sort.by = "none")
  
  # Add miRNA ID column
  res$ID <- fit$genes[match(rownames(res), rownames(logCPMFiltered))]
  res <- res[, c("ID", setdiff(names(res), "ID"))]  # Move ID to the first column
  
  # Add medians
  if ("Group1" %in% strsplit(comp, "vs")[[1]]) {
    res$Median_Group1 <- medianGroup1[selected][match(res$ID, rownames(logCPMFiltered))]
  }
  if ("Group2" %in% strsplit(comp, "vs")[[1]]) {
    res$Median_Group2 <- medianGroup2[selected][match(res$ID, rownames(logCPMFiltered))]
  }
  if ("Group0" %in% strsplit(comp, "vs")[[1]]) {
    res$Median_Group0 <- medianGroup0[selected][match(res$ID, rownames(logCPMFiltered))]
  }

  out_path <- file.path(output_dir, paste0("limma_DEGs_", comp, ".csv"))
  write.table(res, out_path, sep = ",", quote = FALSE, row.names = FALSE)
  message(" Saved: ", out_path)
  }
}

main()

