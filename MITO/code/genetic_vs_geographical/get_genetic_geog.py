#!/usr/bin/env python3
# Works on Python 3.8+
# Merges genetic vs geographic distances for ALL species with robust label cleaning.

import os
import re
from typing import Optional
import pandas as pd

GENETIC_DIR = "genetic_distance"
GEO_DIR = "geographic_distance_matrices"
OUT_DIR = "genetic_vs_geographic"

# Handle multiple naming styles
GENETIC_PATTERNS = [
    "average_genetic_distance_4fold_{sp}.csv",
    "average_sample_distances_{sp}.csv",
    "average_genetic_distance_{sp}.csv",
]

def cut_prefix(s: str, prefix: str) -> str:
    return s[len(prefix):] if s.startswith(prefix) else s

def cut_suffix(s: str, suffix: str) -> str:
    return s[:-len(suffix)] if s.endswith(suffix) else s

def normalize_label(s: str) -> str:
    """Normalize sample labels:
    - strip whitespace and quotes
    - remove a trailing '__suffix' (e.g., '__sexaurelia', '__multimicronucleatum')
    """
    s = str(s).strip().replace('"', '')
    s = re.sub(r'__[^/\\\s]+$', '', s)
    return s

def clean_square_matrix(path: str) -> pd.DataFrame:
    """Read square matrix CSV, clean labels, and align rows/cols to common names."""
    df = pd.read_csv(path, index_col=0)
    df.index = [normalize_label(x) for x in df.index]
    df.columns = [normalize_label(x) for x in df.columns]
    names = sorted(set(df.index) & set(df.columns))
    if not names:
        return pd.DataFrame()
    return df.loc[names, names]

def find_genetic_path(species: str) -> Optional[str]:
    for pat in GENETIC_PATTERNS:
        p = os.path.join(GENETIC_DIR, pat.format(sp=species))
        if os.path.exists(p):
            return p
    return None

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Species from geo filenames
    geo_species = set()
    for fn in os.listdir(GEO_DIR):
        if fn.startswith("geographic_distance_matrix_") and fn.endswith(".csv"):
            sp = cut_suffix(cut_prefix(fn, "geographic_distance_matrix_"), ".csv")
            geo_species.add(sp)

    # Species from genetic filenames (support several prefixes)
    gen_species = set()
    for fn in os.listdir(GENETIC_DIR):
        if not fn.endswith(".csv"):
            continue
        name = cut_suffix(fn, ".csv")
        for prefix in ("average_genetic_distance_4fold_",
                       "average_sample_distances_",
                       "average_genetic_distance_"):
            if name.startswith(prefix):
                gen_species.add(cut_prefix(name, prefix))
                break

    species_list = sorted(geo_species & gen_species)
    if not species_list:
        print("[error] No overlapping species found between the two folders.")
        return

    for sp in species_list:
        geo_path = os.path.join(GEO_DIR, "geographic_distance_matrix_%s.csv" % sp)
        gen_path = find_genetic_path(sp)

        if not (os.path.exists(geo_path) and gen_path and os.path.exists(gen_path)):
            print("[skip] Missing file(s) for %s" % sp)
            continue

        geo = clean_square_matrix(geo_path)
        gen = clean_square_matrix(gen_path)

        if geo.empty or gen.empty:
            print("[warn] %s: empty matrix after cleaning (no common labels)." % sp)
            # Optional quick peek to debug headers
            try:
                gh = list(pd.read_csv(geo_path, nrows=0).columns)[:5]
                ph = list(pd.read_csv(gen_path, nrows=0).columns)[:5]
                print("  geo header peek:", gh)
                print("  gen header peek:", ph)
            except Exception:
                pass
            continue

        samples = sorted(set(geo.index) & set(gen.index))
        if len(samples) < 2:
            print("[warn] %s: <2 overlapping samples; skipping." % sp)
            continue

        geo = geo.loc[samples, samples]
        gen = gen.loc[samples, samples]

        # Long format
        geo_long = (geo.reset_index()
                      .melt(id_vars='index', var_name='Sample2', value_name='Geographic_Distance')
                      .rename(columns={'index': 'Sample1'}))
        gen_long = (gen.reset_index()
                      .melt(id_vars='index', var_name='Sample2', value_name='Genetic_Distance')
                      .rename(columns={'index': 'Sample1'}))

        merged = pd.merge(gen_long, geo_long, on=['Sample1', 'Sample2'], how='inner')

        # Remove self-pairs and collapse (A,B) vs (B,A)
        merged = merged[merged['Sample1'] != merged['Sample2']].copy()
        merged['pair'] = merged.apply(lambda r: tuple(sorted((r['Sample1'], r['Sample2']))), axis=1)
        merged = merged.drop_duplicates(subset='pair').drop(columns='pair')

        out_path = os.path.join(OUT_DIR, "%s_genetic_vs_geographic.csv" % sp)
        merged.to_csv(out_path, index=False)
        print("[ok] %s: %d samples → %d pairs → %s" % (sp, len(samples), len(merged), out_path))

if __name__ == "__main__":
    main()
