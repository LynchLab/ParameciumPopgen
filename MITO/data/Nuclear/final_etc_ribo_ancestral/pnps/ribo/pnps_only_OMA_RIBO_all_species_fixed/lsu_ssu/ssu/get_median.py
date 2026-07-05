import pandas as pd
import glob
import numpy as np

# -----------------------------
# INPUTS
# -----------------------------
csv_pattern = "*.csv"   # run inside ssu folder

# -----------------------------
# combine all CSVs
# -----------------------------
all_data = []

for f in glob.glob(csv_pattern):
    # skip unwanted files if present
    if f.endswith("mean_ssu.csv") or f.endswith(".py"):
        continue

    species = f.replace(".csv", "")

    df = pd.read_csv(f)
    df["Species"] = species

    # use existing OMA_ID as gene
    df["Gene"] = df["OMA_ID"].astype(str)
    df["gene_name"] = df["OMA_ID"].astype(str)

    all_data.append(df)

final = pd.concat(all_data, ignore_index=True)

# -----------------------------
# reorder columns
# -----------------------------
cols = ["Species", "Gene", "OMA_ID", "gene", "pN/pS", "subunit"]
cols = [c for c in cols if c in final.columns] + \
       [c for c in final.columns if c not in cols]

final = final[cols]

# save combined file
final.to_csv("ALL_species_SSU_pnps_with_OMA.csv", index=False)
print("Saved: ALL_species_SSU_pnps_with_OMA.csv")

# -----------------------------
# CLEAN
# -----------------------------
final["pN/pS"] = pd.to_numeric(final["pN/pS"], errors="coerce")
final = final.replace([np.inf, -np.inf], np.nan)
final = final.dropna(subset=["Gene", "pN/pS"])

# keep consistent with mito pipeline
final = final[final["pN/pS"] > 0]

# -----------------------------
# STEP 1: median per gene
# -----------------------------
gene_median = (
    final.groupby("Gene", as_index=False)
         .agg(median_pN_pS=("pN/pS", "median"))
)

# -----------------------------
# STEP 2: overall median
# -----------------------------
overall_median = gene_median["median_pN_pS"].median()

print("SSU nuclear median pN/pS:", round(overall_median, 4))

# optional save
gene_median.to_csv("SSU_nuclear_gene_median_pnps.csv", index=False)
print("Saved: SSU_nuclear_gene_median_pnps.csv")
