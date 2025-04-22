# ===========================
# Import Required Libraries
# ===========================
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted

# ===========================
# Define Input & Output Paths
# ===========================
limma_files = [
    "./reports/tables/limma_DEGs_Group1vsGroup2_filtered.csv",
    "./reports/tables/limma_DEGs_Group2vsGroup0_filtered.csv",
    "./reports/tables/limma_DEGs_Group1vsGroup0_filtered.csv"
]

coefficients_file = "./reports/tables/merged_miRNA_coefficients.csv"
rf_features_file = "./reports/tables/selected_features_per_fold.txt"

output_table = "./reports/tables/combined_DEA_LASSO_RF.csv"
output_venn = "./reports/figures/venn_DEA_LASSO_RF.png"

# ===========================
# Step 1: Load and Combine Unique DEA miRNAs
# ===========================
print(" Loading DEA (LIMMA) files...")
all_mirnas = []

for path in limma_files:
    df = pd.read_csv(path)
    mirna_column = df.iloc[:, 0].dropna()
    all_mirnas.extend(mirna_column.tolist())
    print(f"  ✔ {path}: {len(mirna_column)} miRNAs")

unique_mirnas = sorted(set(all_mirnas))
dea_df = pd.DataFrame(unique_mirnas, columns=["DEA"])

# ===========================
# Step 2: Load LASSO Features
# ===========================
print("\n Loading LASSO feature list...")
coeff_df = pd.read_csv(coefficients_file)
coeff_df["Feature"] = coeff_df["Feature"].str.replace("`", "", regex=False)
lasso_df = pd.DataFrame(coeff_df["Feature"])
lasso_df.columns = ["LASSO"]

# ===========================
# Step 3: Load RF Final Features
# ===========================
print("\n Parsing Random Forest (RF) final features...")
rf_features = []

with open(rf_features_file, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "Final combined features:" in line:
        if i + 1 < len(lines):
            next_line = lines[i + 1]
            rf_features = [feat.strip() for feat in next_line.split(",") if feat.strip()]
        break

rf_df = pd.DataFrame(rf_features, columns=["RF"])
print(f"  ✔ RF features parsed: {len(rf_features)}")

# ===========================
# Step 4: Combine All Columns
# ===========================
max_len = max(len(dea_df), len(lasso_df), len(rf_df))

dea_df = dea_df.reindex(range(max_len)).fillna("")
lasso_df = lasso_df.reindex(range(max_len)).fillna("")
rf_df = rf_df.reindex(range(max_len)).fillna("")

combined_df = pd.concat([dea_df, lasso_df, rf_df], axis=1)

# ===========================
# Step 5: Save Combined Table
# ===========================
combined_df.to_csv(output_table, index=False)
print(f"\n Combined DEA + LASSO + RF saved to:\n{output_table}")

# ===========================
# Step 6: Create Venn Diagram
# ===========================
print("\n Creating Venn diagram...")
set_DEA = set(combined_df["DEA"].dropna().astype(str).str.strip()) - {""}
set_LASSO = set(combined_df["LASSO"].dropna().astype(str).str.strip()) - {""}
set_RF = set(combined_df["RF"].dropna().astype(str).str.strip()) - {""}

plt.figure(figsize=(6, 6))
venn = venn3_unweighted([set_LASSO, set_RF, set_DEA], set_labels=("LASSO", "RF", "DEA"))

plt.title("Venn Diagram of DEA, LASSO, and RF Features (Equal Circle Size)")
plt.tight_layout()
plt.savefig(output_venn, dpi=300)
plt.close()

print(f" Venn diagram saved to: {output_venn}")

# ===========================
# Step 7: Print Shared miRNAs
# ===========================
common_features = set_DEA & set_LASSO & set_RF

print("\n🧬 miRNAs common to DEA, LASSO, and RF:")
if common_features:
    for miRNA in sorted(common_features):
        print(f" - {miRNA}")
else:
    print(" No common miRNAs found.")

print(f"\n Total shared features: {len(common_features)}")

