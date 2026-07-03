#!/usr/bin/env python3

import os

import sys

import argparse

from collections import defaultdict



import pandas as pd

from Bio import SeqIO





def fail(msg: str):

    print(f"ERROR: {msg}", file=sys.stderr)

    sys.exit(1)





def load_unannot_seq_dict(base_dir, species, sample, cache):

    """

    Load and cache unannotated nucleotide sequences for (species, sample).



    Looks for:

        base_dir/species/sample/protein_unannotated/unannotated_nuc.fna



    Returns dict: orf_id -> SeqRecord

    """

    key = (species, sample)

    if key in cache:

        return cache[key]



    fasta_path = os.path.join(

        base_dir, species, sample, "protein_unannotated", "unannotated_nuc.fna"

    )

    if not os.path.isfile(fasta_path):

        print(

            f"[WARN] Missing unannotated_nuc.fna for {species}/{sample}: {fasta_path}",

            file=sys.stderr,

        )

        cache[key] = {}

        return cache[key]



    recs = {r.id: r for r in SeqIO.parse(fasta_path, "fasta")}

    cache[key] = recs

    return recs





def write_cluster_fasta(base_dir, nuc_cache, out_dir, label, members):

    """

    Write FASTA for a cluster.



    label: file label (e.g. 'R1', 'G3', etc.)

    members: list of (species, sample, orf_id)

    """

    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"{label}.fa")



    records = []



    for sp, sa, orf in members:

        seq_dict = load_unannot_seq_dict(base_dir, sp, sa, nuc_cache)

        if orf not in seq_dict:

            print(

                f"[WARN] ORF {orf} not found in unannotated_nuc.fna for {sp}/{sa}",

                file=sys.stderr,

            )

            continue

        rec = seq_dict[orf]

        # informative ID: label|species|sample|orf_id

        rec.id = f"{label}|{sp}|{sa}|{orf}"

        rec.description = ""

        records.append(rec)



    if not records:

        print(

            f"[WARN] No sequences found for cluster {label}; skipping FASTA.",

            file=sys.stderr,

        )

        return



    with open(out_path, "w") as fh:

        SeqIO.write(records, fh, "fasta")



    print(f"[WRITE] {label}: {len(records)} seqs → {out_path}")





def main():

    ap = argparse.ArgumentParser(

        description=(

            "Create dN/dS (global R clusters) and pN/pS (per-species G and R clusters) "

            "FASTA sets from R_clusters.tsv and G_clusters.tsv."

        )

    )

    ap.add_argument(

        "--base-dir",

        required=True,

        help=(

            "Base intergenic_region directory "

            "(contains species/ and RG_clustering/)."

        ),

    )

    args = ap.parse_args()



    base_dir = os.path.abspath(args.base_dir)

    rg_dir = os.path.join(base_dir, "RG_clustering")

    r_clust_tsv = os.path.join(rg_dir, "R_clusters.tsv")

    g_clust_tsv = os.path.join(rg_dir, "G_clusters.tsv")



    if not os.path.isdir(base_dir):

        fail(f"Base dir not found: {base_dir}")

    if not os.path.isdir(rg_dir):

        fail(f"RG_clustering dir not found: {rg_dir}")

    if not os.path.isfile(r_clust_tsv):

        fail(f"Missing R_clusters.tsv: {r_clust_tsv}")

    if not os.path.isfile(g_clust_tsv):

        fail(f"Missing G_clusters.tsv: {g_clust_tsv}")



    print(f"[INFO] Base dir: {base_dir}")

    print(f"[INFO] R clusters file: {r_clust_tsv}")

    print(f"[INFO] G clusters file: {g_clust_tsv}")



    # Read R clusters

    r_df = pd.read_csv(r_clust_tsv, sep="\t", dtype=str)

    required_R = {"species", "sample", "orf_id", "cluster_label"}

    if not required_R.issubset(r_df.columns):

        fail(f"{r_clust_tsv} missing required columns: {required_R - set(r_df.columns)}")



    # Read G clusters

    g_df = pd.read_csv(g_clust_tsv, sep="\t", dtype=str)

    required_G = {"cluster_label", "species", "samples", "orf_ids"}

    if not required_G.issubset(g_df.columns):

        fail(f"{g_clust_tsv} missing required columns: {required_G - set(g_df.columns)}")



    # Root output

    unan_root = os.path.join(base_dir, "unannotated_analysis")

    dnds_root = os.path.join(unan_root, "dnds")

    pnps_root = os.path.join(unan_root, "pnps")



    os.makedirs(unan_root, exist_ok=True)

    os.makedirs(dnds_root, exist_ok=True)

    os.makedirs(pnps_root, exist_ok=True)



    # Cache for nucleotide sequences per (species, sample)

    nuc_cache = {}



    # =========================================================

    # PART 1: dN/dS — global R clusters (multi-species)

    # =========================================================

    print("\n[PART 1] dN/dS: writing global R clusters from R_clusters.tsv...")



    dnds_R_dir = os.path.join(dnds_root, "Rcluster")

    os.makedirs(dnds_R_dir, exist_ok=True)



    if r_df.empty:

        print("[INFO] No rows in R_clusters.tsv.")

    else:

        for clabel, sub in sorted(r_df.groupby("cluster_label"), key=lambda x: x[0]):

            members = [(row["species"], row["sample"], row["orf_id"])

                       for _, row in sub.iterrows()]

            if not members:

                continue

            # Keep global label for dN/dS (e.g. R1, R2, ...)

            write_cluster_fasta(base_dir, nuc_cache, dnds_R_dir, clabel, members)



    # =========================================================

    # PART 2: pN/pS — per-species G and within-species R

    # =========================================================

    print("\n[PART 2] pN/pS: writing per-species G/R clusters...")



    species_from_R = set(r_df["species"].dropna().unique())

    species_from_G = set(g_df["species"].dropna().unique())

    species_list = sorted(species_from_R.union(species_from_G))



    print(f"[INFO] Species found: {', '.join(species_list)}")



    for species in species_list:

        print(f"\n[INFO] Processing species: {species}")



        sp_root = os.path.join(pnps_root, species)

        g_dir = os.path.join(sp_root, "Gcluster")

        r_dir = os.path.join(sp_root, "Rcluster")



        os.makedirs(g_dir, exist_ok=True)

        os.makedirs(r_dir, exist_ok=True)



        # ---------- G clusters (from G_clusters.tsv) ----------

        g_sp = g_df[g_df["species"] == species].copy()

        g_written = 0



        for _, row in g_sp.iterrows():

            label = row["cluster_label"]

            samples_str = row["samples"]

            orfs_str = row["orf_ids"]



            if pd.isna(samples_str) or pd.isna(orfs_str):

                print(

                    f"[WARN] Missing samples/orf_ids for G cluster {label} in {species}",

                    file=sys.stderr,

                )

                continue



            samples = [s.strip() for s in samples_str.split(",") if s.strip()]

            orfs = [o.strip() for o in orfs_str.split(",") if o.strip()]



            if len(samples) != len(orfs):

                print(

                    f"[WARN] G cluster {label} in {species}: "

                    f"len(samples)={len(samples)} != len(orf_ids)={len(orfs)}; "

                    f"using min length.",

                    file=sys.stderr,

                )



            n = min(len(samples), len(orfs))

            if n < 2:

                # for pN/pS you usually want at least 2 seqs

                continue



            members = [(species, samples[i], orfs[i]) for i in range(n)]

            write_cluster_fasta(base_dir, nuc_cache, g_dir, label, members)

            g_written += 1



        # ---------- R clusters within species (from R_clusters.tsv) ----------

        r_sp = r_df[r_df["species"] == species].copy()

        r_written = 0



        for clabel, sub in sorted(r_sp.groupby("cluster_label"), key=lambda x: x[0]):

            samples_set = set(sub["sample"])

            if len(samples_set) < 2:

                # only one sample for this species in this R cluster -> skip for pN/pS

                continue



            members = [(row["species"], row["sample"], row["orf_id"])

                       for _, row in sub.iterrows()]

            if not members:

                continue



            # Use the same R cluster label (e.g. R1) within species

            write_cluster_fasta(base_dir, nuc_cache, r_dir, clabel, members)

            r_written += 1



        print(

            f"[SUMMARY] {species}: "

            f"G clusters written (>=2 seqs) = {g_written}, "

            f"R clusters written (>=2 samples) = {r_written}"

        )



    print("\n[DONE]")

    print(f"  dN/dS R clusters → {dnds_R_dir}")

    print(f"  pN/pS clusters   → {pnps_root}")





if __name__ == "__main__":

    main()


