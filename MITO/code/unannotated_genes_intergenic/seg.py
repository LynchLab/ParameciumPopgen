#!/usr/bin/env python3

import pandas as pd

from collections import defaultdict

import re



CLUST = "unannotated_RGZ_clusters.tsv"

IDMAP = "id_map.tsv"



# Read both as pure strings to avoid weird object types

df = pd.read_csv(CLUST, sep="\t", dtype=str)

idmap = pd.read_csv(IDMAP, sep="\t", dtype=str)



# Build meta: short_id -> (species, sample, orf_id)

meta = {

    row["short_id"]: (row["species"], row["sample"], row["orf_id"])

    for _, row in idmap.iterrows()

}



clean_rows = []



def extract_short_id(text):

    """Extract U0000123-style ID from any corrupted string."""

    m = re.search(r"(U\d{7})", str(text))

    return m.group(1) if m else None



for _, row in df.iterrows():

    sp = str(row["species"])

    sa = str(row["sample"])

    orf = str(row["orf_id"])

    cl = str(row["cluster_label"])



    if "bound method" in sp or "bound method" in sa or "bound method" in orf:

        short = extract_short_id(sp) or extract_short_id(sa) or extract_short_id(orf)

        if short and short in meta:

            sp, sa, orf = meta[short]



    clean_rows.append([str(sp), str(sa), str(orf), str(cl)])



clean_df = pd.DataFrame(clean_rows, columns=["species", "sample", "orf_id", "cluster_label"])



def get_type(label: str) -> str:

    if label.startswith("R"):

        return "R"

    if label.startswith("G"):

        return "G"

    if label.startswith("Z"):

        return "Z"

    return "OTHER"



clean_df["type"] = clean_df["cluster_label"].apply(get_type)



# --------------------------

# R clusters: one row per ORF

# --------------------------

df_R = clean_df[clean_df["type"] == "R"].copy()

df_R = df_R.sort_values(["cluster_label", "species", "sample"], key=lambda col: col.astype(str))

df_R.to_csv("R_clusters.tsv", sep="\t", index=False)



# --------------------------

# Z clusters: one row per ORF

# --------------------------

df_Z = clean_df[clean_df["type"] == "Z"].copy()

df_Z = df_Z.sort_values(["species", "sample", "cluster_label"], key=lambda col: col.astype(str))

df_Z.to_csv("Z_clusters.tsv", sep="\t", index=False)



# --------------------------

# G clusters: one row per cluster, all samples in one field

# --------------------------

G_rows = []

for cl, sub in clean_df[clean_df["type"] == "G"].groupby("cluster_label"):

    species = sorted(sub["species"].unique())[0]  # within-species cluster

    samples = ",".join(sorted(sub["sample"].unique()))

    orfs = ",".join(sorted(sub["orf_id"].unique()))

    G_rows.append([cl, species, samples, orfs])



G_df = pd.DataFrame(G_rows, columns=["cluster_label", "species", "samples", "orf_ids"])

G_df = G_df.sort_values("cluster_label", key=lambda col: col.astype(str))

G_df.to_csv("G_clusters.tsv", sep="\t", index=False)



print("[DONE] Cleaned and wrote:")

print("  R_clusters.tsv")

print("  G_clusters.tsv")

print("  Z_clusters.tsv")


