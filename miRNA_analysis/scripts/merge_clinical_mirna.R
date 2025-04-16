#!/usr/bin/env Rscript

# --------- Load libraries ---------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# --------- Parse arguments ---------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript merge_clinical_mirna.R <clinical_data.csv> <mirna_data.tsv> <output_file.csv>")
}

clinical_path <- args[1]
mirna_path <- args[2]
output_path <- args[3]

# --------- Load data ---------
cat("Loading clinical data...\n")
clinical_data <- read_csv(clinical_path) %>%
  rename(status = Group)

cat("Loading miRNA data...\n")
miRNA_data <- read_tsv(mirna_path) %>%
  as.data.frame()

# --------- Transform miRNA matrix ---------
rownames(miRNA_data) <- miRNA_data$miRNA
miRNA_data <- miRNA_data[, -1]
miRNA_data_transposed <- as.data.frame(t(miRNA_data)) %>%
  rownames_to_column(var = "miRNA")

# --------- Extract and clean IDs ---------
miRNA_data_transposed <- miRNA_data_transposed %>%
  mutate(core_id = str_to_upper(str_extract(miRNA, "[^_]+$")))

clinical_data <- clinical_data %>%
  mutate(core_id = str_to_upper(miRNA))

# Remove non-printable characters
clinical_data$core_id <- str_replace_all(clinical_data$core_id, "[[:cntrl:]]", "")

# Apply manual fixes
manual_fixes <- c(
  "KGKONTROLE16" = "KGKONTR16",
  "KGKONTROLE17" = "KGKONTR17",
  "KGKONTROLE23" = "KGKONTR23",
  "KGKONTROLE39" = "KGKONTR39",
  "KGKONTROLE42" = "KGKONTR42",
  "22.11.K2" = "22.11.K2"
)

clinical_data$core_id <- recode(clinical_data$core_id, !!!manual_fixes)

# --------- Merge ---------
cat(" Merging datasets...\n")
merged_data <- full_join(clinical_data, miRNA_data_transposed, by = "core_id") %>%
  select(-miRNA.y) %>%
  rename(miRNA = miRNA.x)

# --------- Diagnostics ---------
missing_from_mirna <- setdiff(clinical_data$core_id, miRNA_data_transposed$core_id)
missing_from_clinical <- setdiff(miRNA_data_transposed$core_id, clinical_data$core_id)

cat(" Missing from miRNA data:\n")
print(missing_from_mirna)

cat("\n Missing from clinical data:\n")
print(missing_from_clinical)

# --------- Save output ---------
cat("Saving merged data to:", output_path, "\n")
write_csv(merged_data, output_path)

cat("Done!\n")

