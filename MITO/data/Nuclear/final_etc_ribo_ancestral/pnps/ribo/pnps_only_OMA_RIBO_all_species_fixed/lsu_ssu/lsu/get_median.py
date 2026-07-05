import pandas as pd
import glob
import numpy as np

# -----------------------------
# INPUTS
# -----------------------------
csv_pattern = "*.csv"   # run inside lsu folder

# -----------------------------
# combine all CSVs
# -----------------------------
all_data = []

for f in glob.glob(csv_pattern):
    # skip summary/script files if present
    if f in ["nuclear_pnps_mean_lsu.csv", "oma_id.py"]:
        continue

    species = f.replace(".csv", "")

    df = pd.read_csv(f)
    df["Species"] = species

    # use existing OMA_ID as gene name
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
final.to_csv("ALL_species_LSU_pnps_with_OMA.csv", index=False)
print("Saved: ALL_species_LSU_pnps_with_OMA.csv")

# -----------------------------
# compute median pN/pS
# -----------------------------
final["pN/pS"] = pd.to_numeric(final["pN/pS"], errors="coerce")
final = final.replace([np.inf, -np.inf], np.nan)
final = final.dropna(subset=["Gene", "pN/pS"])

# keep same as mito pipeline
final = final[final["pN/pS"] > 0]

# median per OMA gene across species
gene_median = (
    final.groupby("Gene", as_index=False)
         .agg(median_pN_pS=("pN/pS", "median"))
)

# overall median across genes
overall_median = gene_median["median_pN_pS"].median()

print("LSU nuclear median pN/pS:", round(overall_median, 4))

# optional save
gene_median.to_csv("LSU_nuclear_gene_median_pnps.csv", index=False)
print("Saved: LSU_nuclear_gene_median_pnps.csv")
