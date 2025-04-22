import argparse
import pandas as pd
import numpy as np
import os
import pyreadstat
from clinical_utils import (
    report_missing_categorical,
    report_missing_float,
    compute_missing_rates,
    summarize_categorical_stats,
    summarize_numerical_stats,
    inject_missing_values,
    impute_with_median,
    impute_with_iterative,
    impute_with_miceforest,
    impute_with_em,
    impute_with_hotdeck,
    evaluate_imputation,
    plot_error_comparison,
    analyze_dataset,
    scalers, 
    apply_scaler,
    plot_histograms,
    yeo_johnson_transform
)
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")


# ----------------------------- Constants ----------------------------- #
NUMERICAL_COLS = [
    'Age', 'BMI', 'WBC', 'RBC', 'HGB', 'HCT', 'PLT', 'NE_abs', 'LY_abs', 'MO_abs', 'Eo_Abs',
    'Ba_abs', 'Total_cholesterol', 'LDL', 'HDL', 'TG', 'Creatinine', 'ALAT', 'ASAT', 'Glucose',
    'HbA1c', 'Target_fibrotic', 'Target_lipidic', 'Target_necrotic', 'Target_Calcific'
]

CATEGORICAL_COLS = [
    'miRNA', 'Group', 'Sex', 'Smoking', 'Positive_family_history',
    'Art_hipert', 'Previous_PCI', 'Previous_MI'
]


COLUMNS_TO_IMPUTE = [
    'BMI', 'HCT', 'NE_abs', 'LY_abs', 'MO_abs', 'Eo_Abs', 'Ba_abs', 
    'Total_cholesterol', 'HDL', 'TG', 'Creatinine', 'ALAT', 'ASAT', 'Glucose'
]

EXCLUDE_FROM_IMPUTATION = [
    'Target_fibrotic', 'Target_lipidic', 'Target_necrotic', 'Target_Calcific', 'HbA1c'
]

COLUMN_NAME_MAPPING = {
    'NE_abs': 'NE abs', 'LY_abs': 'LY abs', 'MO_abs': 'MO abs', 'Eo_Abs': 'Eo Abs', 'Ba_abs': 'Ba abs',
    'Total_cholesterol': 'Total Cholesterol', 'Target_fibrotic': 'Fibrotic tissue',
    'Target_lipidic': 'Lipidic tissue', 'Target_necrotic': 'Necrotic tissue', 'Target_Calcific': 'Calcific tissue'
}

METHOD_MAPPING = {
    'Median': 'Median Imputation',
    'Iterative Imputer': 'MICE (Iterative Imputation)',
    'miceforest': 'MICE (Random Forest)',
    'EM': 'Expectation-Maximization (EM)',
    'Hot-Deck': 'Hot-Deck Imputation'
}

# -------------------------- Main Functions --------------------------- #
def load_and_preprocess_data(path1, path2):
    data1, _ = pyreadstat.read_sav(path1)
    data2, _ = pyreadstat.read_sav(path2)
    
    cat1 = [c for c in CATEGORICAL_COLS if c in data1.columns]
    cat2 = [c for c in CATEGORICAL_COLS if c in data2.columns]
    num1 = [c for c in NUMERICAL_COLS if c in data1.columns]
    num2 = [c for c in NUMERICAL_COLS if c in data2.columns]

    data1_filtered = data1[cat1 + num1].copy()
    data2_filtered = data2[cat2 + num2].copy()
    data2_filtered['Group'] = 0

    for df, cats in zip([data1_filtered, data2_filtered], [cat1, cat2]):
        for col in cats:
            df[col] = df[col].astype('category')

    data1_filtered['Positive_family_history'].replace(['', 'unknown', 'missing', 'NA'], pd.NA, inplace=True)
    return data1_filtered, data2_filtered

def clean_and_merge(data1, data2):
    data1 = data1.dropna(subset=['Group'])
    data2 = data2.dropna(subset=['Group'])
    
    data1 = data1[~data1['miRNA'].isin(['KG46', 'KG33', 'KP059'])]
    data2 = data2[~data2['miRNA'].isin(['11K3', '11K1', '9K1', '27K1', '30K1', '31K3', '31K5', '1K3', '2K3', ''])]
    
    if (dup := data2[data2['miRNA'] == '12K1']).shape[0] > 1:
        drop_idx = dup.isnull().sum(axis=1).idxmax()
        data2 = data2.drop(index=drop_idx)

    data2['miRNA'].replace({'20K3': '10.10.Plume', '8K1': '26K1'}, inplace=True)
    return pd.concat([data1, data2], ignore_index=True)

def fix_known_issues(df):
    corrections = {
        (97, 11): 5.2, (117, 11): 3.1, (115, 16): 1.46, (91, 26): np.nan,
        (72, 30): 10.17, (129, 9): 35.7, (136, 12): 155, (101, 6): 0.0
    }
    for (row, col), val in corrections.items():
        df.iat[row, col] = val
    return df

def summarize_and_analyze(df):
    print(f"\n📊 Dataset shape: {df.shape}")
    print(report_missing_categorical(df))
    print(report_missing_float(df))
    print(compute_missing_rates(df, exclude_columns=['miRNA']))
    print("\n🔁 Duplicated rows:", df.duplicated().sum())
    print(summarize_categorical_stats(df))
    print(summarize_numerical_stats(df))

def run_imputation_workflow(df):
    df_missing, missing_idx, orig_missing = inject_missing_values(df, NUMERICAL_COLS, seed=42)
    methods = {
        'Median': impute_with_median(df_missing, NUMERICAL_COLS),
        'Iterative Imputer': impute_with_iterative(df_missing, NUMERICAL_COLS),
        'miceforest': impute_with_miceforest(df_missing, NUMERICAL_COLS),
        'EM': impute_with_em(df_missing, NUMERICAL_COLS),
        'Hot-Deck': impute_with_hotdeck(df_missing, NUMERICAL_COLS)
    }
    mae, rmse = {}, {}
    for method, df_imp in methods.items():
        mae[method], rmse[method] = evaluate_imputation(orig_missing, df_imp, missing_idx)
    mae_df, rmse_df = pd.DataFrame(mae), pd.DataFrame(rmse)
    return mae_df, rmse_df, methods

def select_best_method(mae_df, rmse_df):
    combined = mae_df.mean() + rmse_df.mean()
    return combined.idxmin()

def run_normalization_pipeline(df, best_method):
    numeric_cols = [col for col in df.select_dtypes(include=[np.number]).columns if col != 'Group']
    all_res, normal_feats = analyze_dataset(df, numeric_cols, 'Imputed Dataset', scalers)
    best_norm = None
    if not normal_feats.empty and 'Normal Methods' in normal_feats.columns:
        best_norm = normal_feats['Normal Methods'].str.split(', ').explode().value_counts().idxmax()
    return best_norm, normal_feats


def final_transform(df, best_imputer_fn, norm_method):
    df_imputed = best_imputer_fn(df, COLUMNS_TO_IMPUTE, exclude=EXCLUDE_FROM_IMPUTATION)
    numeric_cols = [col for col in df_imputed.select_dtypes(include=[np.number]).columns if col != 'Group']

    if norm_method and norm_method in scalers:
        for col in numeric_cols:
            try:
                df_imputed[col] = apply_scaler(scalers[norm_method], df_imputed[col])
            except Exception:
                continue

    for col in numeric_cols:
        try:
            df_imputed[col] = yeo_johnson_transform(df_imputed[col])
        except Exception:
            continue

    return df_imputed


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True)
    parser.add_argument('--path2', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default='./reports')
    parser.add_argument('--data_save_path', type=str, default='./data/processed')
    return parser.parse_args()

def main():
    args = parse_args()

    figures_path = os.path.join(args.output_dir, 'figures')
    os.makedirs(figures_path, exist_ok=True)

    df1, df2 = load_and_preprocess_data(args.path1, args.path2)
    df_combined = clean_and_merge(df1, df2)
    df_combined = fix_known_issues(df_combined)

    for col in CATEGORICAL_COLS:
        if col in df_combined.columns:
            df_combined[col] = df_combined[col].astype('category')

    for col in df_combined.select_dtypes(include='category').columns:
        df_combined[col] = df_combined[col].cat.remove_unused_categories()

    # Save cleaned df_combined
    os.makedirs(args.data_save_path, exist_ok=True)
    df_combined.to_csv(os.path.join(args.data_save_path, 'df_combined_cleaned.csv'), index=False)

    # Save original summaries
    tables_path = os.path.join(args.output_dir, 'tables')
    os.makedirs(tables_path, exist_ok=True)
    summarize_categorical_stats(df_combined).to_csv(os.path.join(tables_path, 'categorical_summary_before_imputation.csv'))
    summarize_numerical_stats(df_combined).to_csv(os.path.join(tables_path, 'numerical_summary_before_imputation.csv'))

    summarize_and_analyze(df_combined)

    print("\n🚧 Running Imputation Comparison...\n")
    mae_df, rmse_df, all_methods = run_imputation_workflow(df_combined)

    print("📊 Mean Absolute Error (MAE):")
    print(mae_df.mean().round(4).sort_values())

    print("\n📊 Root Mean Squared Error (RMSE):")
    print(rmse_df.mean().round(4).sort_values())

    best_method = select_best_method(mae_df, rmse_df)
    print(f"\n✅ Best imputation method based on MAE + RMSE: {best_method}")

    best_imputer_fn = {
        'Median': impute_with_median,
        'Iterative Imputer': impute_with_iterative,
        'miceforest': impute_with_miceforest,
        'EM': impute_with_em,
        'Hot-Deck': impute_with_hotdeck
    }[best_method]

    df_imputed = best_imputer_fn(df_combined, COLUMNS_TO_IMPUTE, exclude=EXCLUDE_FROM_IMPUTATION)
    print(f"\n✅ Imputation completed on original df_combined using the best method.")
    print("🔎 Remaining missing values:", df_imputed[NUMERICAL_COLS].isna().sum().sum())

    summarize_and_analyze(df_imputed)

    # Save summaries and imputed data
    df_imputed.to_csv(os.path.join(args.data_save_path, 'df_imputed.csv'), index=False)
    summarize_categorical_stats(df_imputed).to_csv(os.path.join(tables_path, 'categorical_summary_after_imputation.csv'))
    summarize_numerical_stats(df_imputed).to_csv(os.path.join(tables_path, 'numerical_summary_after_imputation.csv'))

    print("\n🔎 Running Normalization Comparison...")
    best_norm, _ = run_normalization_pipeline(df_combined, best_method)

    if best_norm:
        print(f"\n🏆 Best normalization method based on normality: {best_norm}")
    else:
        print("⚠️ No transformations produced normal distributions for any variable.")

    df_final = final_transform(df_combined.copy(), best_imputer_fn, best_norm)
    print("✅ Final normalization applied to imputed dataset.")

    # Save normalized data
    df_final.to_csv(os.path.join(args.data_save_path, 'df_normalized.csv'), index=False)

    plot_error_comparison(
        mae_df.rename(index=COLUMN_NAME_MAPPING).rename(columns=METHOD_MAPPING),
        'MAE Comparison',
        'MAE',
        'mae_comparison.png',
        output_dir=figures_path
    )

    plot_error_comparison(
        rmse_df.rename(index=COLUMN_NAME_MAPPING).rename(columns=METHOD_MAPPING),
        'RMSE Comparison',
        'RMSE',
        'rmse_comparison.png',
        output_dir=figures_path
    )

    plot_histograms(df_combined, all_methods[best_method], df_final, output_path=os.path.join(figures_path, 'histogram_numerical_original_imputed_transformed.png'))
    print(f"📁 Plot saved to: {os.path.join(figures_path, 'histogram_numerical_original_imputed_transformed.png')}")

    mae_df.to_csv(os.path.join(tables_path, 'mae_values.csv'))
    rmse_df.to_csv(os.path.join(tables_path, 'rmse_values.csv'))

    print("📁 Tables saved to:", tables_path)
    print("✅ Pipeline completed.")

if __name__ == '__main__':
    main()

