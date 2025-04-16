# ===========================
# 0. Argument Parsing
# ===========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript filter_limma_results.R <output_dir> <input_file1> [<input_file2> ...]")
}

output_dir <- args[1]
input_files <- args[-1]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ===========================
# 1. Filtering Function
# ===========================
filter_limma_results <- function(input_path, output_dir, logFC_threshold = log2(2), pval_threshold = 0.05) {
  data <- read.csv(input_path, sep = ",")
  
  # Filter based on logFC and adjusted p-value
  filtered <- data[abs(data$logFC) > logFC_threshold & data$adj.P.Val < pval_threshold, ]
  
  # Create output file name
  input_filename <- basename(input_path)
  output_filename <- sub("\\.csv$", "_filtered.csv", input_filename)
  output_path <- file.path(output_dir, output_filename)
  
  # Save filtered results
  write.csv(filtered, output_path, row.names = FALSE)
  
  cat("Filtered file saved to:", output_path, "\n")
}

# ===========================
# 2. Process Each Input File
# ===========================
for (file in input_files) {
  filter_limma_results(file, output_dir)
}

