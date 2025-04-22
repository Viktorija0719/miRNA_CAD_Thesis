# ===========================
# 0. Setup and Argument Parsing
# ===========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript miRNA_RF_feature_selection.R <input_csv> <figures_dir> <tables_dir>")
}

input_file <- args[1]
figures_dir <- args[2]
tables_dir <- args[3]

# Create output directories if they don't exist
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

# ===========================
# 1. Load Libraries
# ===========================
library(readr)
library(dplyr)
library(caret)
library(randomForest)
library(ggplot2)

set.seed(123)

# ===========================
# 2. Load and Prepare Data
# ===========================
df <- read_csv(input_file)
miRNA_cols <- names(df)[grepl("^hsa", names(df))]
df$status <- as.factor(df$status)

train_idx <- createDataPartition(df$status, p = 0.7, list = FALSE)
df_train <- df[train_idx, ]
df_test  <- df[-train_idx, ]

x_train <- df_train[, miRNA_cols]
x_test  <- df_test[, miRNA_cols]
y_train <- df_train$status
y_test  <- df_test$status

# ===========================
# 3. Feature Selection (5-Fold CV)
# ===========================
k <- 5
folds <- createFolds(y_train, k = k, list = TRUE)

feature_selection_results <- data.frame(Fold = integer(), Num_Features = integer(), OOB_Error = numeric())
all_selected_features <- list()
importance_aggregated <- data.frame(Feature = character(), MeanDecreaseGini = numeric())

for (i in seq_along(folds)) {
  cat("Processing fold", i, "\n")
  
  test_idx <- folds[[i]]
  train_fold <- df_train[-test_idx, ]
  y_fold <- train_fold$status
  x_fold <- train_fold[, miRNA_cols]
  
  rf_initial <- randomForest(x = x_fold, y = y_fold, ntree = 750, importance = TRUE)
  
  importance_df <- as.data.frame(importance(rf_initial))
  importance_df$Feature <- rownames(importance_df)
  importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
  importance_aggregated <- rbind(importance_aggregated, importance_df[, c("Feature", "MeanDecreaseGini")])
  
  top_features <- importance_df$Feature
  current_features <- top_features[1:5]
  best_oob_error <- tail(rf_initial$err.rate[, "OOB"], 1)
  
  for (f in top_features[-(1:5)]) {
    temp_features <- c(current_features, f)
    temp_model <- randomForest(x = train_fold[, temp_features], y = y_fold, ntree = 750, importance = TRUE)
    temp_oob <- tail(temp_model$err.rate[, "OOB"], 1)
    
    if (temp_oob < best_oob_error) {
      best_oob_error <- temp_oob
      current_features <- temp_features
    }
  }
  
  feature_selection_results <- rbind(feature_selection_results, data.frame(
    Fold = i,
    Num_Features = length(current_features),
    OOB_Error = best_oob_error
  ))
  all_selected_features[[paste0("Fold_", i)]] <- current_features
}

# ===========================
# 4. Final Model Training and Evaluation
# ===========================
all_features_combined <- unique(unlist(all_selected_features))

rf_final <- randomForest(x = df_train[, all_features_combined], y = y_train, ntree = 750, importance = TRUE)
rf_predictions <- predict(rf_final, newdata = df_test[, all_features_combined])
final_conf_matrix <- confusionMatrix(rf_predictions, y_test)

# ===========================
# 5. Save Outputs
# ===========================

# Save selected features per fold
features_file <- file.path(tables_dir, "selected_features_per_fold.txt")
sink(features_file)
for (i in seq_along(all_selected_features)) {
  cat(paste0("Fold ", i, ": "), paste(all_selected_features[[i]], collapse = ", "), "\n")
}
cat("\nFinal combined features:\n")
cat(paste(all_features_combined, collapse = ", "), "\n")
sink()

# Save confusion matrix
conf_matrix_file <- file.path(tables_dir, "confusion_matrix.txt")
sink(conf_matrix_file)
print(final_conf_matrix)
sink()

# Save feature importance summary table
importance_summary <- importance_aggregated %>%
  group_by(Feature) %>%
  summarise(MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(MeanDecreaseGini))

final_importance <- importance_summary %>%
  filter(Feature %in% all_features_combined)

write_csv(final_importance, file.path(tables_dir, "final_feature_importance.csv"))

# Save bar plot of feature importance
plot_file <- file.path(figures_dir, "feature_importance_plot.png")
png(plot_file, width = 900, height = 600)
ggplot(final_importance, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top miRNA Features by Mean Decrease in Gini",
    x = "miRNA Feature",
    y = "Mean Decrease in Gini"
  ) +
  theme_minimal()
dev.off()

