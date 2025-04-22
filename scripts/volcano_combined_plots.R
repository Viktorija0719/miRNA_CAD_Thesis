
#!/usr/bin/env Rscript

# Load required libraries
library(EnhancedVolcano)
library(ggplot2)
library(patchwork)
library(cowplot)

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 4) {
    stop("Usage: Rscript volcano_combined_plots.R <Group1vsGroup2.csv> <Group2vsGroup0.csv> <Group1vsGroup0.csv> <output_plot.png>")
  }
  
  file1 <- argv[1]
  file2 <- argv[2]
  file3 <- argv[3]
  output_file <- argv[4]
  
  # Load data
  res1 <- read.table(file1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  res2 <- read.table(file2, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  res3 <- read.table(file3, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # Top 10 miRNAs per comparison
  top_miRNAs1 <- res1[order(res1$adj.P.Val), ]$ID[1:10]
  top_miRNAs2 <- res2[order(res2$adj.P.Val), ]$ID[1:10]
  top_miRNAs3 <- res3[order(res3$adj.P.Val), ]$ID[1:10]
  
  # Volcano plot creator
  create_volcano <- function(data, title, top_miRNAs) {
    EnhancedVolcano(data,
                    lab = ifelse(data$ID %in% top_miRNAs, data$ID, ""),
                    x = 'logFC',
                    y = 'adj.P.Val',
                    xlim = c(-6, 5),
                    ylim = c(0, 15),
                    title = title,
                    subtitle = NULL,
                    caption = NULL,
                    xlab = "Log2 Fold Change",
                    ylab = "-Log10 P-value",
                    pCutoff = 0.05,
                    FCcutoff = log2(2),
                    pointSize = 3.0,
                    labSize = 5
    ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none"
      )
  }
  
  # Create volcano plots
  p1 <- create_volcano(res1, "A", top_miRNAs1)
  p2 <- create_volcano(res2, "B", top_miRNAs2)
  p3 <- create_volcano(res3, "C", top_miRNAs3)
  
  # Extract legend from a dummy plot
  legend_plot <- EnhancedVolcano(res1,
                                 lab = ifelse(res1$ID %in% top_miRNAs1, res1$ID, ""),
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                                 xlim = c(-6, 5),
                                 ylim = c(0, 15),
                                 pCutoff = 0.05,
                                 FCcutoff = log2(2),
                                 pointSize = 3.0,
                                 labSize = 5
  ) +
    theme(legend.position = "right")
  
  legend <- cowplot::get_legend(legend_plot)
  
  # Arrange layout
  final_plot <- (p1 + p2) / (p3 + plot_grid(legend, NULL, ncol = 1, rel_heights = c(1, 0))) +
    plot_layout(heights = c(1, 1))
  
  # Save final plot
  ggsave(output_file, plot = final_plot, width = 12, height = 10, dpi = 300)
}

main()

