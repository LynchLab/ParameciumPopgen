#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(".")
LSU_DIR = BASE / "lsu"
SSU_DIR = BASE / "ssu"

# =========================================================
# CLEAN FUNCTION
# =========================================================
def clean_series(s):
    s = pd.to_numeric(s, errors="coerce")
    s = s.replace([np.inf, -np.inf], np.nan)
    return s.dropna()

# =========================================================
# 1. OVERALL LSU + SSU
# =========================================================
all_lsu = []
all_ssu = []

for f in LSU_DIR.glob("*.csv"):
    df = pd.read_csv(f)
    if "pN/pS" not in df.columns:
        continue
    all_lsu.extend(clean_series(df["pN/pS"]))

for f in SSU_DIR.glob("*.csv"):
    df = pd.read_csv(f)
    if "pN/pS" not in df.columns:
        continue
    all_ssu.extend(clean_series(df["pN/pS"]))

all_lsu = pd.Series(all_lsu)
all_ssu = pd.Series(all_ssu)

print("\n================ OVERALL =================")
print(f"LSU  -> mean: {all_lsu.mean():.4f}, median: {all_lsu.median():.4f}, n={len(all_lsu)}")
print(f"SSU  -> mean: {all_ssu.mean():.4f}, median: {all_ssu.median():.4f}, n={len(all_ssu)}")

# =========================================================
# 2. PER SPECIES
# =========================================================
results = []

species_files = sorted(set([f.name for f in LSU_DIR.glob("*.csv")]))

for fname in species_files:
    lsu_file = LSU_DIR / fname
    ssu_file = SSU_DIR / fname

    species = fname.replace(".csv", "")

    # LSU
    if lsu_file.exists():
        df_lsu = pd.read_csv(lsu_file)
        lsu_vals = clean_series(df_lsu["pN/pS"]) if "pN/pS" in df_lsu.columns else pd.Series()
    else:
        lsu_vals = pd.Series()

    # SSU
    if ssu_file.exists():
        df_ssu = pd.read_csv(ssu_file)
        ssu_vals = clean_series(df_ssu["pN/pS"]) if "pN/pS" in df_ssu.columns else pd.Series()
    else:
        ssu_vals = pd.Series()

    results.append({
        "species": species,
        "LSU_mean": lsu_vals.mean(),
        "LSU_median": lsu_vals.median(),
        "LSU_n": len(lsu_vals),
        "SSU_mean": ssu_vals.mean(),
        "SSU_median": ssu_vals.median(),
        "SSU_n": len(ssu_vals),
    })

res_df = pd.DataFrame(results)

print("\n================ PER SPECIES =================")
print(res_df.to_string(index=False, float_format="%.4f"))

# save
res_df.to_csv("LSU_SSU_per_species_summary.csv", index=False)

print("\nSaved: LSU_SSU_per_species_summary.csv")
