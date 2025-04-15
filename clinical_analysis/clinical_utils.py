import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, shapiro, ttest_ind, mannwhitneyu
import numpy as np
import miceforest as mf
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.mixture import GaussianMixture
import scipy.spatial.distance
from sklearn.metrics import mean_absolute_error, mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, MaxAbsScaler, PowerTransformer, QuantileTransformer
from scipy.stats import shapiro, anderson, yeojohnson
from statsmodels.stats.diagnostic import lilliefors
from scipy.stats import mstats
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def report_missing_categorical(df):
    """Report missing values for categorical columns."""
    cat_missing = df.select_dtypes(include='category').isna().sum()
    missing_cols = cat_missing[cat_missing > 0]

    if not missing_cols.empty:
        total_rows = len(df)
        percent_missing = (missing_cols / total_rows) * 100
        summary_df = pd.DataFrame({
            'Missing Count': missing_cols,
            'Percent Missing': percent_missing
        }).sort_values(by='Percent Missing', ascending=False)
        return summary_df
    return pd.DataFrame(columns=['Missing Count', 'Percent Missing'])


def report_missing_float(df):
    """Report missing values for float columns."""
    float_missing = df.select_dtypes(include='float64').isna().sum()
    missing_cols = float_missing[float_missing > 0]

    if not missing_cols.empty:
        total_rows = len(df)
        percent_missing = (missing_cols / total_rows) * 100
        summary_df = pd.DataFrame({
            'Missing Count': missing_cols,
            'Percent Missing': percent_missing
        }).sort_values(by='Percent Missing', ascending=False)
        return summary_df
    return pd.DataFrame(columns=['Missing Count', 'Percent Missing'])


def compute_missing_rates(df, exclude_columns=None):
    """
    Compute total, row-wise, and column-wise missing rates.
    Excludes columns like 'miRNA' if specified.
    """
    if exclude_columns:
        df = df.drop(columns=exclude_columns, errors='ignore')

    missing_matrix = df.isnull().astype(int)
    total = missing_matrix.sum().sum()
    total_missing_rate = total / (missing_matrix.shape[0] * missing_matrix.shape[1])

    row_missing_rate = (
        missing_matrix.sum(axis=1).apply(lambda x: 1 if x > 0 else 0).sum()
        / missing_matrix.shape[0]
    )

    column_missing_rate = (
        missing_matrix.sum(axis=0).apply(lambda x: 1 if x > 0 else 0).sum()
        / missing_matrix.shape[1]
    )

    return {
        "Total Missing Rate": total_missing_rate,
        "Row Missing Rate": row_missing_rate,
        "Column Missing Rate": column_missing_rate
    }



def summarize_categorical_stats(df, group_col='Group'):
    """Generate counts, percentages, and pairwise p-values for categorical variables."""
    df = df.copy()
    df['Positive_family_history'] = pd.to_numeric(df['Positive_family_history'], errors='coerce')

    categorical_columns = ['Sex', 'Smoking', 'Positive_family_history', 'Art_hipert', 'Previous_PCI', 'Previous_MI']
    group_pairs = [(1.0, 0.0), (2.0, 0.0), (1.0, 2.0)]

    all_result_keys = []
    for col in categorical_columns:
        if col == 'Smoking':
            all_result_keys.extend([f'{col}_1 (Smoking)', f'{col}_2 (Quitted)'])
        else:
            all_result_keys.append(col)

    results = pd.DataFrame(index=all_result_keys)

    # Count + percent summary per group
    for col in categorical_columns:
        if col not in df.columns:
            continue

        if col == 'Smoking':
            for val, label in zip([1, 2], ['Smoking', 'Quitted']):
                key = f'{col}_{val} ({label})'
                for group in sorted(df[group_col].dropna().unique()):
                    group_subset = df[df[group_col] == group]
                    count = (group_subset[col] == val).sum()
                    percent = (group_subset[col] == val).mean() * 100
                    results.at[key, f'Group {group}'] = f"{count} ({percent:.1f}%)"
        else:
            for group in sorted(df[group_col].dropna().unique()):
                group_subset = df[df[group_col] == group]
                counts = group_subset[col].value_counts()
                percents = group_subset[col].value_counts(normalize=True) * 100
                formatted = [
                    f"{int(count)} ({percent:.1f}%)"
                    for count, percent in zip(counts, percents)
                ]
                results.at[col, f'Group {group}'] = ' | '.join(formatted)

    # Pairwise p-values
    for col in categorical_columns:
        for g1, g2 in group_pairs:
            subset = df[df[group_col].isin([g1, g2])]
            contingency = pd.crosstab(subset[group_col], subset[col])

            try:
                if contingency.shape == (2, 2):
                    _, p_val = fisher_exact(contingency)
                else:
                    _, p_val, _, _ = chi2_contingency(contingency)
            except Exception:
                p_val = float('nan')

            label = f'p-Value Group {g1} vs Group {g2}'

            if col == 'Smoking':
                results.at['Smoking_1 (Smoking)', label] = f"{p_val:.2e}" if pd.notna(p_val) else "NaN"
                results.at['Smoking_2 (Quitted)', label] = f"{p_val:.2e}" if pd.notna(p_val) else "NaN"
            else:
                results.at[col, label] = f"{p_val:.2e}" if pd.notna(p_val) else "NaN"

    return results





def summarize_numerical_stats(df, group_col='Group'):
    """Generate descriptive stats and perform pairwise comparisons between groups for numerical variables."""
    numerical_columns = [
        'Age', 'BMI', 'WBC', 'RBC', 'HGB', 'HCT', 'PLT',
        'NE_abs', 'LY_abs', 'MO_abs', 'Eo_Abs', 'Ba_abs', 'Total_cholesterol',
        'LDL', 'HDL', 'TG', 'Creatinine', 'ALAT', 'ASAT', 'Glucose',
        'HbA1c', 'Target_fibrotic', 'Target_lipidic', 'Target_necrotic', 'Target_Calcific'
    ]

    df[group_col] = pd.to_numeric(df[group_col], errors='coerce')
    group_pairs = [(1.0, 0.0), (2.0, 0.0), (1.0, 2.0)]

    results_df = pd.DataFrame(index=numerical_columns, columns=[
        'Group 0.0', 'Group 1.0', 'Group 2.0',
        'Stat Group 1.0 vs Group 0.0', 'p-Value Group 1.0 vs Group 0.0',
        'Stat Group 2.0 vs Group 0.0', 'p-Value Group 2.0 vs Group 0.0',
        'Stat Group 1.0 vs Group 2.0', 'p-Value Group 1.0 vs Group 2.0'
    ])

    for group, group_data in df.groupby(group_col):
        for col in numerical_columns:
            data = group_data[col].dropna()
            if len(data) > 3:
                stat, p = shapiro(data)
                if p > 0.05:
                    results_df.loc[col, f'Group {group}'] = f"{data.mean():.2f} ± {data.std():.2f}"
                else:
                    results_df.loc[col, f'Group {group}'] = f"{data.median():.2f} (IQR: {(data.quantile(0.75) - data.quantile(0.25)):.2f})"
            else:
                results_df.loc[col, f'Group {group}'] = "NaN"

    for col in numerical_columns:
        for g1, g2 in group_pairs:
            d1 = df[df[group_col] == g1][col].dropna()
            d2 = df[df[group_col] == g2][col].dropna()
            if len(d1) > 3 and len(d2) > 3:
                p1 = shapiro(d1)[1]
                p2 = shapiro(d2)[1]
                if p1 > 0.05 and p2 > 0.05:
                    stat, p_val = ttest_ind(d1, d2)
                    stat_str = f"t={stat:.2f}"
                else:
                    stat, p_val = mannwhitneyu(d1, d2)
                    stat_str = f"U={stat:.2f}"
                if (g1, g2) == (1.0, 0.0):
                    results_df.loc[col, 'Stat Group 1.0 vs Group 0.0'] = stat_str
                    results_df.loc[col, 'p-Value Group 1.0 vs Group 0.0'] = f"{p_val:.2e}"
                elif (g1, g2) == (2.0, 0.0):
                    results_df.loc[col, 'Stat Group 2.0 vs Group 0.0'] = stat_str
                    results_df.loc[col, 'p-Value Group 2.0 vs Group 0.0'] = f"{p_val:.2e}"
                elif (g1, g2) == (1.0, 2.0):
                    results_df.loc[col, 'Stat Group 1.0 vs Group 2.0'] = stat_str
                    results_df.loc[col, 'p-Value Group 1.0 vs Group 2.0'] = f"{p_val:.2e}"
            else:
                label_map = {
                    (1.0, 0.0): 'Group 1.0 vs Group 0.0',
                    (2.0, 0.0): 'Group 2.0 vs Group 0.0',
                    (1.0, 2.0): 'Group 1.0 vs Group 2.0'
                }
                label = label_map[(g1, g2)]
                results_df.loc[col, f'Stat {label}'] = "Insufficient data"
                results_df.loc[col, f'p-Value {label}'] = "Insufficient data"

    return results_df







# ------------------------- Missing Value Injection ------------------------- #
def inject_missing_values(df, columns, n_missing=10, seed=42):
    np.random.seed(seed)
    df_copy = df.copy()
    missing_indices = {}
    original_missing_values = {}

    for column in columns:
        non_missing = df_copy[df_copy[column].notna()].index
        selected = np.random.choice(non_missing, size=n_missing, replace=False)
        missing_indices[column] = selected
        original_missing_values[column] = df_copy.loc[selected, column].copy()
        df_copy.loc[selected, column] = np.nan

    return df_copy, missing_indices, original_missing_values


# ------------------------- Imputation Methods ------------------------- #
def impute_with_median(df, columns, exclude=None):
    df_copy = df.copy()
    if exclude:
        columns = [col for col in columns if col not in exclude]

    for col in columns:
        median = df_copy[col].median()
        df_copy[col] = df_copy[col].fillna(median)
    return df_copy


def impute_with_iterative(df, columns, exclude=None):
    df_copy = df.copy()
    if exclude:
        columns = [col for col in columns if col not in exclude]

    imputer = IterativeImputer(max_iter=10, random_state=10, sample_posterior=True)
    df_copy[columns] = imputer.fit_transform(df_copy[columns])
    return df_copy


def impute_with_miceforest(df, columns, exclude=None):
    df_copy = df.copy()
    if exclude:
        columns = [col for col in columns if col not in exclude]

    kernel = mf.ImputationKernel(df_copy[columns], datasets=1, random_state=10, save_all_iterations=True)
    kernel.mice(max_iterations=10)
    df_imputed = kernel.complete_data(dataset=0)
    df_copy[columns] = df_imputed[columns]
    return df_copy


def impute_with_em(df, columns, exclude=None):
    df_copy = df.copy()
    if exclude:
        columns = [col for col in columns if col not in exclude]

    for col in columns:
        missing = df_copy[col].isnull()
        if missing.sum() > 0:
            gm = GaussianMixture(n_components=1, random_state=10)
            known = df_copy.loc[~missing, col].values.reshape(-1, 1)
            gm.fit(known)
            imputed = gm.sample(missing.sum())[0].flatten()
            df_copy.loc[missing, col] = imputed
    return df_copy



def impute_with_hotdeck(df, columns, exclude=None):
    df_copy = df.copy()
    if exclude:
        columns = [col for col in columns if col not in exclude]

    numeric = df_copy.select_dtypes(include=[np.number])
    rows_with_nan = numeric[numeric[columns].isnull().any(axis=1)].copy()
    rows_without_nan = numeric.dropna()

    for i in range(len(rows_with_nan)):
        row = rows_with_nan.iloc[i, :]
        null_mask = row.isnull()
        row_clean = row.loc[~null_mask].values.reshape(1, -1)
        donor_matrix = rows_without_nan.loc[:, ~null_mask].values
        distances = scipy.spatial.distance.cdist(row_clean, donor_matrix, 'euclidean')
        donor_idx = np.argmin(distances)

        for col in columns:
            if null_mask[col]:
                value = rows_without_nan.iloc[donor_idx, numeric.columns.get_loc(col)]
                rows_with_nan.iloc[i, numeric.columns.get_loc(col)] = value

    df_copy.update(rows_with_nan)
    return df_copy


# ------------------------- Evaluation ------------------------- #
def evaluate_imputation(original, imputed, missing_indices):
    mae_scores = {}
    rmse_scores = {}

    for column, indices in missing_indices.items():
        orig_vals = original[column]
        imputed_vals = imputed.loc[indices, column]

        valid = orig_vals.notna() & imputed_vals.notna()
        orig_vals = orig_vals[valid]
        imputed_vals = imputed_vals[valid]

        if len(orig_vals) > 0:
            mae_scores[column] = mean_absolute_error(orig_vals, imputed_vals)
            rmse_scores[column] = np.sqrt(mean_squared_error(orig_vals, imputed_vals))

    return mae_scores, rmse_scores


# ------------------------- Plotting ------------------------- #
def plot_error_comparison(df, title, ylabel, filename, output_dir="E:/RSU_work/Plots", custom_palette=None):  #!!!!!!
    df_melted = df.reset_index().melt(id_vars='index', var_name='Imputation Method', value_name=ylabel)

    plt.figure(figsize=(15, 8))
    sns.barplot(data=df_melted, x='index', y=ylabel, hue='Imputation Method', 
                palette=custom_palette if custom_palette else 'viridis')

    plt.xlabel('Columns', fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title='Imputation Method', fontsize=12)

    # Save the plot
    filepath = f"{output_dir}/PLT_{filename}"
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to: {filepath}")










# ------------------------- Normalization + Normality ------------------------- #
# Define additional transformation methods
def min_max_normalize(series):
    return (series - series.min()) / (series.max() - series.min())

def max_abs_scaling(series):
    return series / series.abs().max()

def tanh_normalization(series):
    return np.tanh(series)

def z_score_normalize(series):
    return (series - series.mean()) / series.std()

def robust_scale(series):
    median = series.median()
    IQR = series.quantile(0.75) - series.quantile(0.25)
    return (series - median) / IQR

def quantile_transform(data):
    sorted_data = np.sort(data)
    ranks = np.argsort(np.argsort(data))
    quantiles = (ranks + 1) / (len(data) + 1)
    return quantiles

def square_root_transform(series):
    return np.sqrt(series)

def square_transformation(series):
    return np.square(series)

def exponential_transformation(series, exponent=0.5):
    return np.power(series, exponent)

def yeo_johnson_transform(series):
    transformed_data, _ = yeojohnson(series)
    return transformed_data

def l1_normalization(series):
    return series / np.sum(np.abs(series))

def l2_normalization(series):
    return series / np.sqrt(np.sum(np.square(series)))

def decimal_scaling(series):
    j = np.ceil(np.log10(series.abs().max()))
    return series / (10**j)

def winsorize_transformation(series, limits=[0.05, 0.05]):
    return mstats.winsorize(series, limits=limits)

# Define transformation methods including the new ones
scalers = {
    'Original': None,
    'Logarithmic': lambda x: np.log1p(x),   
    'Min-Max Normalization': min_max_normalize,
    'Min-Max Scaling': MinMaxScaler(),
    'Z-Score Normalization': z_score_normalize,
    'Robust Scaling': robust_scale,
    'Max-Abs Scaling': max_abs_scaling, 
    'Power Transformation(Yeo-Johnson)': yeo_johnson_transform,   
    'Quantile Transformation': QuantileTransformer(output_distribution='normal'),
    'Winsorized': lambda x: mstats.winsorize(x, limits=[0.05, 0.05]),
    'Tanh Normalization': tanh_normalization,
    'Square Root Transform': square_root_transform,
    'Square Transformation': square_transformation,
    'Exponential Transformation': exponential_transformation,
    'L1 Normalization': l1_normalization,
    'L2 Normalization': l2_normalization,
    'Decimal Scaling': decimal_scaling,
}



def apply_scaler(scaler, data):
    if scaler is None:
        return data
    elif isinstance(scaler, (MinMaxScaler, StandardScaler, RobustScaler, MaxAbsScaler, PowerTransformer, QuantileTransformer)):
        return scaler.fit_transform(data.values.reshape(-1, 1)).flatten()
    else:
        return scaler(data)


def check_normality(data):
    p_values = {}
    try:
        lillie_stat, lillie_p = lilliefors(data, dist='norm')
        shapiro_stat, shapiro_p = shapiro(data)
        anderson_stat = anderson(data)
        anderson_result = 1 if anderson_stat.statistic < anderson_stat.critical_values[2] else 0

        p_values['Lilliefors'] = lillie_p
        p_values['Shapiro'] = shapiro_p
        p_values['Anderson'] = anderson_result

    except Exception:
        p_values = {'Lilliefors': np.nan, 'Shapiro': np.nan, 'Anderson': np.nan}

    return p_values


def analyze_dataset(data, columns_to_analyze, dataset_name, scalers_dict):
    all_results = []
    normal_methods = {}

    for column in columns_to_analyze:
        normal_methods[column] = []
        for method_name, scaler in scalers_dict.items():
            try:
                transformed = apply_scaler(scaler, data[column].dropna())

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    warnings.simplefilter("ignore", category=UserWarning)

                    lillie_p = lilliefors(transformed, dist='norm')[1]
                    anderson_test = anderson(transformed)
                    anderson_passed = int(anderson_test.statistic < anderson_test.critical_values[2])
                    shapiro_p = shapiro(transformed)[1]

                passed = sum([
                    lillie_p > 0.05 if not np.isnan(lillie_p) else False,
                    shapiro_p > 0.05 if not np.isnan(shapiro_p) else False,
                    anderson_passed
                ])

                all_results.append({
                    'Dataset': dataset_name,
                    'Column': column,
                    'Transformation': method_name,
                    'Lilliefors': lillie_p,
                    'Anderson': 'Normal' if anderson_passed else 'Not Normal',
                    'Shapiro': shapiro_p,
                    'Normality Assumed': 'Yes' if passed >= 2 else 'No'
                })

                if passed >= 2:
                    normal_methods[column].append(method_name)

            except Exception as e:
                continue

    # Final formatting
    all_results_df = pd.DataFrame(all_results)
    normal_features_df = pd.DataFrame([
        (col, ', '.join(methods)) for col, methods in normal_methods.items() if methods
    ], columns=['Column', 'Normal Methods'])

    return all_results_df, normal_features_df


def plot_histograms(original, df_imputed, df_transformed, output_path=None):

    def check_normality(data):
        p_values = []
        p_values.append(lilliefors(data, dist='norm')[1])
        anderson_test = anderson(data)
        p_values.append(anderson_test.statistic < anderson_test.critical_values[2])  # Anderson-Darling
        p_values.append(shapiro(data)[1])
        return p_values

    numeric_columns = original.select_dtypes(include=[np.number]).columns
    plt.figure(figsize=(18, len(numeric_columns) * 3))
    pastel_colors = {'green': '#77dd77', 'blue': '#75bbfd'}

    for i, column in enumerate(numeric_columns, start=1):
        original_data = original[column].dropna()
        imputed_data = df_imputed[column].dropna()
        transformed_data = df_transformed[column].dropna()

        def assign_color(p_vals):
            return pastel_colors['green'] if sum(p > 0.05 for p in [p_vals[0], p_vals[2]]) + int(p_vals[1]) >= 2 else pastel_colors['blue']

        # Plot original
        plt.subplot(len(numeric_columns), 3, 3 * i - 2)
        plt.hist(original_data, bins=20, color=assign_color(check_normality(original_data)), edgecolor='black')
        plt.title(f'{column} (Original)')
        plt.xlabel('Value')
        plt.ylabel('Count')

        # Plot imputed
        plt.subplot(len(numeric_columns), 3, 3 * i - 1)
        plt.hist(imputed_data, bins=20, color=assign_color(check_normality(imputed_data)), edgecolor='black')
        plt.title(f'{column} (Imputed)')
        plt.xlabel('Value')
        plt.ylabel('Count')

        # Plot transformed
        plt.subplot(len(numeric_columns), 3, 3 * i)
        plt.hist(transformed_data, bins=20, color=assign_color(check_normality(transformed_data)), edgecolor='black')
        plt.title(f'{column} (Yeo-Johnson)')
        plt.xlabel('Value')
        plt.ylabel('Count')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, format='png', bbox_inches='tight')
        print(f"📁 Plot saved to: {output_path}")
    #plt.show()

