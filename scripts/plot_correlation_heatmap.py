import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib.colors import LinearSegmentedColormap

def parse_args():
    parser = argparse.ArgumentParser(description="Plot correlation heatmap from CSV results.")
    parser.add_argument("--input", type=str, required=True, help="Path to the input correlation CSV file")
    parser.add_argument("--output", type=str, required=True, help="Path to save the output heatmap image")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load data
    df = pd.read_csv(args.input)

    # Rename mapping
    rename_map = {
        'Age': 'Age', 'Sex': 'Sex', 'BMI': 'BMI', 'WBC': 'WBC', 'RBC': 'RBC', 'HGB': 'HGB', 'HCT': 'HCT', 'PLT': 'PLT',
        'NE_abs': 'NE abs', 'LY_abs': 'LY abs', 'MO_abs': 'MO abs', 'Eo_Abs': 'Eo Abs', 'Total_cholesterol': 'Total Cholesterol',
        'LDL': 'LDL', 'HDL': 'HDL', 'TG': 'TG', 'Creatinine': 'Creatinine', 'ALAT': 'ALAT', 'ASAT': 'ASAT',
        'Glucose': 'Glucose', 'HbA1c': 'HbA1c', 'Target_fibrotic': 'Fibrotic tissue', 'Target_lipidic': 'Lipidic tissue',
        'Target_necrotic': 'Necrotic tissue', 'Target_Calcific': 'Calcific tissue',
        'status': 'Status', 'Smoking_1': 'Smoking', 'Smoking_2': 'Quitted smoking',
        'Positive_family_history': 'Positive family history', 'Art_hipert': 'Arterial Hypertension',
        'Previous_PCI': 'Previous PCI', 'Previous_MI': 'Previous MI'
    }

    # Rename clinical variables
    df['Variable1'] = df['Variable1'].replace(rename_map)

    # Filter based on thresholds
    filtered_df = df[(abs(df["Correlation Coefficient"]) > 0.20) & (df["p.adjusted"] < 0.05)]

    # Pivot to heatmap structure
    heatmap_data = filtered_df.pivot(index="Variable2", columns="Variable1", values="Correlation Coefficient")
    heatmap_data = heatmap_data.dropna(axis=1, how='all')

    # Retain only desired clinical variable columns in order
    clinical_variables_order = [var for var in rename_map.values() if var in heatmap_data.columns]
    heatmap_data = heatmap_data[clinical_variables_order]

    # Create custom colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_colormap", ["#8e82fe", "white", "#ef4026"])

    # Plot
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        heatmap_data,
        cmap=custom_cmap,
        annot=False,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={'label': 'Correlation Coefficient'},
        square=True
    )
    plt.xlabel("Clinical Variables", fontsize=14)
    plt.ylabel("miRNAs", fontsize=14)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    # Save plot
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f" Heatmap saved to: {args.output}")

    plt.show()

if __name__ == "__main__":
    main()

