#!/usr/bin/env Rscript

# ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(optparse)
})

# ---- Argument Parsing ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input merged data CSV"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./", help = "Output directory [default: %default]"),
  make_option(c("-p", "--prefix"), type = "character", default = "correlation", help = "Prefix for output files [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load and validate data ----
merged_data <- read_csv(opt$input)

hsa_vars <- names(merged_data)[grepl("^hsa-", names(merged_data))]
continuous_vars <- c('Age', 'BMI', 'WBC', 'RBC', 'HGB', 'HCT', 'PLT', 
                     'NE_abs', 'LY_abs', 'MO_abs', 'Eo_Abs', 'Total_cholesterol', 
                     'LDL', 'HDL', 'TG', 'Creatinine', 'ALAT', 'ASAT', 'Glucose', 'HbA1c', 
                     'Target_fibrotic', 'Target_lipidic', 'Target_necrotic', 'Target_Calcific')
discrete_vars <- c('status', 'Sex', 'Smoking_1', 'Smoking_2', 'Positive_family_history', 'Art_hipert', 'Previous_PCI', 'Previous_MI')

# ---- Functions ----
calculate_kendall_correlation_hsa <- function(data, cont_vars, hsa_vars) {
  results <- expand.grid(Variable1 = cont_vars, Variable2 = hsa_vars) %>%
    rowwise() %>%
    mutate(
      test = list(cor.test(data[[Variable1]], data[[Variable2]], method = "kendall")),
      Tau = test$estimate,
      p.value = test$p.value
    ) %>%
    ungroup() %>%
    select(-test) %>%
    mutate(p.adjusted = p.adjust(p.value, method = "BH"))
  return(results)
}

calculate_point_biserial_correlation_hsa <- function(data, hsa_vars, disc_vars) {
  results <- list()
  for (disc in disc_vars) {
    if (length(unique(data[[disc]])) != 2) next
    binary_var <- as.numeric(as.factor(data[[disc]])) - 1
    for (hsa in hsa_vars) {
      test <- cor.test(binary_var, data[[hsa]], method = "pearson")
      results[[paste(disc, hsa, sep = "_")]] <- data.frame(
        Variable1 = disc, Variable2 = hsa,
        Rpb = test$estimate, p.value = test$p.value
      )
    }
  }
  df <- bind_rows(results)
  df$p.adjusted <- p.adjust(df$p.value, method = "BH")
  return(df)
}

combine_results <- function(kendall, biserial) {
  bind_rows(
    kendall %>%
      rename(`Correlation Coefficient` = Tau) %>%
      mutate(Method = "Kendall"),
    
    biserial %>%
      rename(`Correlation Coefficient` = Rpb) %>%
      mutate(Method = "Point-Biserial")
  )
}



# ---- Run analysis ----
cat(" Running Kendall and Point-Biserial correlations...\n")
kendall <- calculate_kendall_correlation_hsa(merged_data, continuous_vars, hsa_vars)
biserial <- calculate_point_biserial_correlation_hsa(merged_data, hsa_vars, discrete_vars)
combined <- combine_results(kendall, biserial)

# ---- Save results ----
output_file <- file.path(opt$output_dir, paste0(opt$prefix, "_combined_results.csv"))
write_csv(combined, output_file)
cat(" Correlation results saved to:", output_file, "\n")

