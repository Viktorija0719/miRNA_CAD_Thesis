# =========================================
# Multinomial LASSO with miRNA (CLI version)
# =========================================

library(glmnet)
library(caret)
library(pROC)
library(readr)
library(dplyr)

# ============================
# Handle Command-Line Arguments
# ============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file.csv> <output_dir/>")
}
input_file <- args[1]
output_dir <- args[2]

# ============================
# Load and Prepare Data
# ============================
set.seed(123)
df <- read_csv(input_file)

# Select miRNA features only
miRNA_cols <- names(df)[grepl("^hsa", names(df))]
X_miRNA <- model.matrix(~ . - 1, df[, miRNA_cols])
y <- as.factor(df$status)

# Outer 5-fold CV
outer_folds <- 5
outer_folds_idx <- createFolds(y, k = outer_folds, list = TRUE)

selected_features_list <- list()
overall_conf_matrices <- list()

for (i in 1:outer_folds) {
  cat("\n========== Fold", i, "==========\n")
  
  test_idx <- outer_folds_idx[[i]]
  train_idx <- setdiff(seq_along(y), test_idx)
  
  X_train <- X_miRNA[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X_miRNA[test_idx, ]
  y_test <- y[test_idx]
  
  alpha_values <- seq(1, 1, by = 0.1)
  best_cvm <- Inf
  best_alpha <- NA
  best_lambda <- NA
  best_model <- NULL
  
  for (alpha_val in alpha_values) {
    set.seed(123)
    lasso_cv <- cv.glmnet(X_train, y_train,
                          family = "multinomial", alpha = alpha_val,
                          type.multinomial = "grouped", nfolds = 5,
                          type.measure = "class")
    
    mean_cvm <- min(lasso_cv$cvm)
    cat("Alpha:", alpha_val, "| Mean Accuracy:", round(1 - mean_cvm, 4), "\n")
    
    if (mean_cvm < best_cvm) {
      best_cvm <- mean_cvm
      best_alpha <- alpha_val
      best_lambda <- lasso_cv$lambda.min
      best_model <- lasso_cv
    }
  }
  
  # Extract selected features
  coefs_list <- coef(best_model, s = best_lambda)
  selected_features <- unique(unlist(lapply(coefs_list, function(coef_matrix) {
    rownames(coef_matrix)[which(coef_matrix != 0)]
  })))
  selected_features <- setdiff(selected_features, "(Intercept)")
  selected_features_list[[i]] <- selected_features
  
  # Fit final model on selected features
  X_train_final <- X_train[, selected_features, drop = FALSE]
  X_test_final <- X_test[, selected_features, drop = FALSE]
  
  final_model <- glmnet(X_train_final, y_train,
                        family = "multinomial",
                        type.multinomial = "grouped",
                        alpha = best_alpha,
                        lambda = best_lambda)
  
  pred_class <- predict(final_model, newx = X_test_final, s = "lambda.min", type = "class")
  pred_class <- factor(pred_class, levels = levels(y_test))
  
  cm <- confusionMatrix(pred_class, y_test)
  print(cm)
  overall_conf_matrices[[paste0("Fold_", i)]] <- cm
}

# ============================
# Performance Summary
# ============================

# Overall accuracy
overall_accuracy <- sapply(overall_conf_matrices, function(cm) cm$overall["Accuracy"])
cat("\n====== Average Accuracy Across Folds ======\n")
print(round(mean(overall_accuracy), 4))

# By-class metrics
extract_by_class_metrics <- function(cm) as.data.frame(cm$byClass)
all_by_class <- lapply(overall_conf_matrices, extract_by_class_metrics)
by_class_combined <- do.call(rbind, all_by_class)

metric_names <- colnames(all_by_class[[1]])
by_class_matrix <- sapply(metric_names, function(metric) {
  sapply(1:3, function(class_idx) {
    mean(sapply(all_by_class, function(m) m[class_idx, metric]), na.rm = TRUE)
  })
})
by_class_matrix <- t(by_class_matrix)
rownames(by_class_matrix) <- metric_names
colnames(by_class_matrix) <- c("Class_0", "Class_1", "Class_2")

cat("\n====== MEAN STATISTICS BY CLASS (ACROSS FOLDS) ======\n")
print(round(by_class_matrix, 4))

# Overall metrics
extract_overall <- function(cm) cm$overall
overall_df <- do.call(rbind, lapply(overall_conf_matrices, extract_overall))
mean_overall <- colMeans(overall_df, na.rm = TRUE)

cat("\n====== MEAN OVERALL STATISTICS (ACROSS FOLDS) ======\n")
print(round(mean_overall, 4))

# ============================
# Extract Coefficients (Final Fold)
# ============================
get_non_zero_coefs_multinomial <- function(model, lambda_value) {
  coefs_list <- lapply(coef(model, s = lambda_value), as.matrix)
  lapply(coefs_list, function(coef_matrix) {
    nz <- coef_matrix[coef_matrix != 0, , drop = FALSE]
    data.frame(Feature = rownames(nz), Coefficient = nz[, 1])
  })
}

non_zero_coefs_final <- get_non_zero_coefs_multinomial(best_model, best_lambda)

# Merge into one table
coef_class_0 <- non_zero_coefs_final[["0"]][, c("Feature", "Coefficient")]
coef_class_1 <- non_zero_coefs_final[["1"]][, c("Feature", "Coefficient")]
coef_class_2 <- non_zero_coefs_final[["2"]][, c("Feature", "Coefficient")]

colnames(coef_class_0) <- c("Feature", "Class_0")
colnames(coef_class_1) <- c("Feature", "Class_1")
colnames(coef_class_2) <- c("Feature", "Class_2")

merged_coefs <- full_join(coef_class_0, coef_class_1, by = "Feature") %>%
  full_join(coef_class_2, by = "Feature") %>%
  filter(Feature != "(Intercept)")

# ============================
# Save Outputs
# ============================
write.csv(merged_coefs, file.path(output_dir, "merged_miRNA_coefficients.csv"), row.names = FALSE)
write.csv(round(by_class_matrix, 4), file.path(output_dir, "metrics_by_class.csv"))
write.csv(round(mean_overall, 4), file.path(output_dir, "metrics_overall.csv"))

cat("\n✅ Files saved to:", output_dir, "\n")

