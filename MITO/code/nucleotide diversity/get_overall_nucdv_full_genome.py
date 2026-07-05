#!/usr/bin/env python3

import os

import sys

import glob

import math

import pandas as pd



# Reads per-site CSVs produced by get_site_nuc_div.py and summarizes:

# - overall pi per species (site-weighted)

# - pi at 1-fold/2-fold/3-fold/4-fold sites

#

# Usage:

#   python3 summarize_species_pi.py /path/to/pn_ps_new

#

# Expects each species directory to contain:

#   per_site_pi/*.csv   (including *_per_site_pi_ALL_GENES.csv OR per-gene per-site CSVs)



OUTNAME = "ALL_species_overall_pi_summary.csv"



def safe_mean(series: pd.Series) -> float:

    s = pd.to_numeric(series, errors="coerce")

    s = s[~s.isna()]

    if len(s) == 0:

        return float("nan")

    return float(s.mean())



def summarize_one_species(species_dir: str) -> dict | None:

    species = os.path.basename(os.path.abspath(species_dir))

    per_site_dir = os.path.join(species_dir, "per_site_pi")

    if not os.path.isdir(per_site_dir):

        return None



    # Prefer the combined ALL_GENES file if present

    combined = glob.glob(os.path.join(per_site_dir, "*_per_site_pi_ALL_GENES.csv"))

    if combined:

        df = pd.read_csv(combined[0])

    else:

        # Otherwise read all per-gene per-site csvs

        files = sorted(glob.glob(os.path.join(per_site_dir, "*_per_site_pi.csv")))

        if not files:

            return None

        df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)



    # Ensure required columns exist

    req = {"pi_site", "is_1fold", "is_2fold", "is_3fold", "is_4fold"}

    missing = req - set(df.columns)

    if missing:

        raise RuntimeError(f"{species}: missing columns in per-site tables: {sorted(missing)}")



    # Overall pi across ALL valid sites (exclude NaN pi_site)

    pi_all = safe_mean(df["pi_site"])

    n_sites_all = int(pd.to_numeric(df["pi_site"], errors="coerce").notna().sum())



    # Per degeneracy class pi (site-weighted within class)

    def class_stats(flag_col: str):

        sub = df[pd.to_numeric(df[flag_col], errors="coerce") == 1].copy()

        pi = safe_mean(sub["pi_site"])

        n = int(pd.to_numeric(sub["pi_site"], errors="coerce").notna().sum())

        return pi, n



    pi_1, n1 = class_stats("is_1fold")

    pi_2, n2 = class_stats("is_2fold")

    pi_3, n3 = class_stats("is_3fold")

    pi_4, n4 = class_stats("is_4fold")



    return {

        "species": species,

        "pi_overall_all_sites": pi_all,

        "n_valid_sites_all": n_sites_all,

        "pi_1fold": pi_1, "n_valid_1fold_sites": n1,

        "pi_2fold": pi_2, "n_valid_2fold_sites": n2,

        "pi_3fold": pi_3, "n_valid_3fold_sites": n3,

        "pi_4fold": pi_4, "n_valid_4fold_sites": n4,

    }



def main():

    root = sys.argv[1] if len(sys.argv) > 1 else "."

    root = os.path.abspath(root)



    # Species dirs = immediate subfolders of root

    species_dirs = sorted([p for p in glob.glob(os.path.join(root, "*")) if os.path.isdir(p)])



    rows = []

    for spdir in species_dirs:

        try:

            r = summarize_one_species(spdir)

        except Exception as e:

            print(f"⚠️  Skipping {spdir}: {e}")

            continue

        if r is None:

            continue

        rows.append(r)



    if not rows:

        sys.exit("❌ No species summaries generated. Did you run get_site_nuc_div.py to create per_site_pi/ ?")



    out_df = pd.DataFrame(rows).sort_values("species")

    out_csv = os.path.join(root, OUTNAME)

    out_df.to_csv(out_csv, index=False)

    print(f"✅ Wrote: {out_csv}")

    print(out_df.to_string(index=False))



if __name__ == "__main__":

    main()


