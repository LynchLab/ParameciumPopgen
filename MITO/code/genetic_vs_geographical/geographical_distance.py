#!/usr/bin/env python3
# Compute per-species geographic distance matrices from a single CSV.
# Expected columns: species, sample, latitude, longitude

import os
import re
import numpy as np
import pandas as pd
from geopy.distance import geodesic

# ==== CONFIG (edit these two lines if you want different paths) ====
INPUT_CSV = "315_geographic_cordinate.csv"
OUTDIR = "geographic_distance_matrices"
# ===================================================================

def safe_name(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(s).strip())

def compute_distance_matrix(df_species: pd.DataFrame) -> pd.DataFrame:
    df_species = df_species.reset_index(drop=True)
    n = len(df_species)
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        ci = (df_species.loc[i, 'latitude'], df_species.loc[i, 'longitude'])
        for j in range(i, n):
            cj = (df_species.loc[j, 'latitude'], df_species.loc[j, 'longitude'])
            d = geodesic(ci, cj).kilometers
            M[i, j] = d
            M[j, i] = d
    return pd.DataFrame(M, index=df_species['sample'], columns=df_species['sample'])

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)
    required = {"species", "sample", "latitude", "longitude"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # Clean types / rows
    df = df.dropna(subset=["species", "sample", "latitude", "longitude"]).copy()
    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    df = df.dropna(subset=["latitude", "longitude"])
    df = df.drop_duplicates(subset=["species", "sample"])

    species_list = sorted(df["species"].unique())

    for sp in species_list:
        sub = df[df["species"] == sp].copy()
        if sub.empty:
            continue
        dist_df = compute_distance_matrix(sub)
        out_name = f"geographic_distance_matrix_{safe_name(sp)}.csv"
        out_path = os.path.join(OUTDIR, out_name)
        dist_df.to_csv(out_path)
        print(f"[ok] {sp}: {len(sub)} samples → {out_path}")

if __name__ == "__main__":
    main()


