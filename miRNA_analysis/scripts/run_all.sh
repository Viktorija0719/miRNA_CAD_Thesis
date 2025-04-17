#!/bin/bash

# ============================================
# Full miRNA and Clinical Data Analysis Script
# ============================================

# NOTE:
# Before running this script, make sure to:
# 1. Complete the Snakemake pipeline for raw miRNA processing
# 2. Preprocess clinical data using:
#    python3 ./clinical_analysis/clinical_main.py  \
#      --path1 "./data/clinical/AAR_control_cases_anon.sav" \
#      --path2 "./data/clinical/AAR_kontroles.sav" \
#      --output_dir "./reports/" \
#      --data_save_path "./data/clinical/"


# ================
# Step 1: Preprocess miRNA counts
# ================
python3 ./miRNA_analysis/scripts/1_mergeMiRPrecursors.py \
  ./data/processed/ \
  ./data/metadata/CeGaT_id_External_id_9.csv \
  ./data/processed/countTable.txt

Rscript ./miRNA_analysis/scripts/2_table_to_matrix.R \
  ./data/processed/countTable.txt \
  ./data/processed/countData.txt

python3 ./miRNA_analysis/scripts/3_filter_out_low_counts.py \
  ./data/processed/countData.txt \
  ./reports/tables/kept_medians.txt \
  ./data/processed/kept_data.txt

Rscript ./miRNA_analysis/scripts/4_combat_seq_correct.R \
  ./data/processed/kept_data.txt \
  ./data/processed/ComBat_corrected_counts.txt

Rscript ./miRNA_analysis/scripts/5_normalize_counts_tmm.R \
  ./data/processed/ComBat_corrected_counts.txt \
  ./data/processed/tmm_normalized_data.csv

# ================
# Step 2: Merge with clinical data
# ================
Rscript ./miRNA_analysis/scripts/merge_clinical_mirna.R \
  ./data/clinical/df_normalized.csv \
  ./data/processed/tmm_normalized_data.csv \
  ./data/processed/merged_clinical_mirna.csv

# ================
# Step 3: Exploratory Analysis
# ================
Rscript ./miRNA_analysis/scripts/plot_pca_miRNA.R \
  ./data/processed/tmm_normalized_data.csv \
  ./reports/figures/pca_plot.png

Rscript ./miRNA_analysis/scripts/correlation_analysis.R \
  --input ./data/processed/merged_clinical_mirna.csv \
  --output_dir ./reports/tables/ \
  --prefix mirna_corr

python3 ./miRNA_analysis/scripts/plot_correlation_heatmap.py \
  --input ./reports/tables/mirna_corr_combined_results.csv \
  --output ./reports/figures/miRNA_heatmap_0.20_0.05.png


# ================
# Step 4: Differential Expression & Visualization
# ================
Rscript ./miRNA_analysis/scripts/limma_DE_analysis.R \
  ./data/processed/tmm_normalized_data.csv \
  ./reports/tables/

Rscript ./miRNA_analysis/scripts/volcano_combined_plots.R \
  ./reports/tables/limma_DEGs_Group1vsGroup2.csv \
  ./reports/tables/limma_DEGs_Group2vsGroup0.csv \
  ./reports/tables/limma_DEGs_Group1vsGroup0.csv \
  ./reports/figures/volcano_combined.png

Rscript ./miRNA_analysis/scripts/filter_limma_results.R \
  ./reports/tables/ \
  ./reports/tables/limma_DEGs_Group1vsGroup2.csv \
  ./reports/tables/limma_DEGs_Group2vsGroup0.csv \
  ./reports/tables/limma_DEGs_Group1vsGroup0.csv


# ================
# Step 5: Feature Selection
# ================
Rscript ./miRNA_analysis/scripts/lasso_multiclass_miRNA.R \
  ./data/processed/merged_clinical_mirna.csv \
  ./reports/tables/

Rscript ./miRNA_analysis/scripts/RF_feature_selection.R \
  ./data/processed/merged_clinical_mirna.csv \
  ./reports/figures/ \
  ./reports/tables/

python3 ./miRNA_analysis/scripts/venn_compare_feature_sets.py


echo "All analysis steps completed successfully!"
