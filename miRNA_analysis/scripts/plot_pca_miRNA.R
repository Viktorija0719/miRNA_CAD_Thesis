#!/usr/bin/env Rscript

# DESCRIPTION:
# Perform PCA on a normalized miRNA count matrix and plot PC1 vs PC2.
# 
# USAGE:
#   Rscript plot_pca_miRNA.R input_file.csv output_plot.png

library(ggplot2)
library(scales)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {

  if (length(argv) != 2) {
    stop("Usage: Rscript plot_pca_miRNA.R <input_file.csv> <output_plot.png>")
  }

  input_file <- argv[1]
  output_file <- argv[2]

  # Load normalized expression matrix (CSV format, comma-separated)
  dat <- read.delim(input_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  # Extract sample names and transpose matrix
  sample_names <- colnames(dat)[-1]
  dat2 <- t(dat[, -1])

  # Perform PCA
  dat.pca <- prcomp(dat2, center = TRUE, scale. = FALSE)
  explained_variance <- summary(dat.pca)$importance[2, 1:2] * 100

  # Assign group labels based on sample name prefix
  sample_groups <- ifelse(grepl("^00_CAD|^0_CAD", sample_names), "Control",
                          ifelse(grepl("^1_CAD", sample_names), "Early CAD",
                                 ifelse(grepl("^2_CAD", sample_names), "Late CAD", "Unknown")))
  sample_groups <- factor(sample_groups, levels = c("Control", "Early CAD", "Late CAD", "Unknown"))

  # Prepare data for plotting
  pca_df <- data.frame(
    Sample = sample_names,
    PC1 = dat.pca$x[, 1],
    PC2 = dat.pca$x[, 2],
    Group = sample_groups
  )

  # Plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Control" = "orange",
    "Early CAD" = "coral3",
    "Late CAD" = "cyan4"
  )) +
  stat_ellipse(type = "t", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black", linetype = "blank"),
    legend.key = element_rect(fill = "white", linetype = "blank"),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  labs(
    title = "",
    x = paste0("PC1 (", sprintf("%.2f", explained_variance[1]), "% variance)"),
    y = paste0("PC2 (", sprintf("%.2f", explained_variance[2]), "% variance)"),
    color = NULL
  )

  # Save the plot
  ggsave(output_file, plot = p, width = 7, height = 5, dpi = 300)
}

main()
