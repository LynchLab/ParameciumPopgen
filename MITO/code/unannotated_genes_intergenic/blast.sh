#!/bin/bash

#SBATCH --job-name=RGZ_blast

#SBATCH --output=RGZ_blast.out

#SBATCH --error=RGZ_blast.err

#SBATCH --time=24:00:00

#SBATCH --cpus-per-task=8

#SBATCH --mem=32G



set -euo pipefail



############################################

# ENVIRONMENT

############################################

module load mamba/latest

module load blast-plus-2.12.0-6z

source activate pyenv



############################################

# CONFIG

############################################



BASE="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/intergenic_region"



WORK="${BASE}/RG_clustering"

EVALUE_MAX="1e-3"

IDENTITY_MIN="30.0"

COV_MIN="0.60"



RAW_FASTA="${WORK}/raw_unannotated.faa"

ALL_FASTA="${WORK}/all_unannotated.faa"

IDMAP_TSV="${WORK}/id_map.tsv"

BLAST_DB="${WORK}/all_unannot_db"

BLAST_OUT="${WORK}/all_unannot_blast.tsv"



mkdir -p "$WORK"



echo "Base dir: $BASE"

echo "Work dir: $WORK"

echo



############################################

# STEP 1 — Collect all unannotated_protein.faa

############################################



echo "[STEP 1] Collecting all unannotated_protein.faa → $RAW_FASTA"

> "$RAW_FASTA"



for species_dir in "$BASE"/*; do

    [ -d "$species_dir" ] || continue

    species=$(basename "$species_dir")



    for sample_dir in "$species_dir"/*; do

        [ -d "$sample_dir" ] || continue

        sample=$(basename "$sample_dir")



        faa="${sample_dir}/protein_unannotated/unannotated_protein.faa"

        [ -s "$faa" ] || continue



        echo "  Adding: $species / $sample"



        awk -v sp="$species" -v sa="$sample" '

            /^>/ { sub(/^>/,""); print ">" sp "|" sa "|" $0; next }

            { print }

        ' "$faa" >> "$RAW_FASTA"

    done

done



############################################

# STEP 2 — Convert headers to short IDs

############################################



echo

echo "[STEP 2] Converting headers to BLAST-safe short IDs..."



export RAW_FASTA ALL_FASTA IDMAP_TSV



python3 << 'EOF'

import os



raw = os.environ["RAW_FASTA"]

allf = os.environ["ALL_FASTA"]

mapf = os.environ["IDMAP_TSV"]



counter = 1



with open(raw) as fin, open(allf, "w") as fout, open(mapf, "w") as fm:

    fm.write("short_id\tspecies\tsample\torf_id\n")



    for line in fin:

        if line.startswith(">"):

            header = line[1:].strip()

            parts = header.split("|", 2)

            if len(parts) == 3:

                species, sample, orf = parts

            else:

                species, sample, orf = "NA", "NA", header



            short = f"U{counter:07d}"

            counter += 1



            fm.write(f"{short}\t{species}\t{sample}\t{orf}\n")

            fout.write(f">{short}\n")

        else:

            fout.write(line)



print("[INFO] BLAST-safe FASTA:", allf)

print("[INFO] ID map:", mapf)

EOF



############################################

# STEP 3 — BLAST DB + ALL-vs-ALL BLASTP

############################################



echo

echo "[STEP 3] Building BLAST DB..."

makeblastdb -in "$ALL_FASTA" -dbtype prot -parse_seqids -out "$BLAST_DB"



echo "[STEP 3] Running BLASTP all-vs-all..."

blastp -query "$ALL_FASTA" -db "$BLAST_DB" -evalue "$EVALUE_MAX" -outfmt "6 qseqid sseqid pident length qlen slen evalue" -num_threads 8 -out "$BLAST_OUT"

echo "[STEP 3] BLAST output: $BLAST_OUT"



############################################

# STEP 4 — Cluster & assign R/G/Z labels

############################################



echo

echo "[STEP 4] Clustering and assigning R/G/Z labels..."



export BLAST_OUT ALL_FASTA IDMAP_TSV IDENTITY_MIN COV_MIN EVALUE_MAX



python3 << 'EOF'

import os

import pandas as pd

from collections import defaultdict



blast = os.environ["BLAST_OUT"]

idmap = os.environ["IDMAP_TSV"]

IDENT = float(os.environ["IDENTITY_MIN"])

COV = float(os.environ["COV_MIN"])

EVAL = float(os.environ["EVALUE_MAX"])



# ---- read id map ----

id_df = pd.read_csv(idmap, sep="\t", dtype=str)

meta = {r.short_id: (r.species, r.sample, r.orf_id) for _, r in id_df.iterrows()}

all_ids = set(meta.keys())



# ---- read BLAST ----

cols = ["q", "s", "pident", "alen", "qlen", "slen", "evalue"]

try:

    df = pd.read_csv(blast, sep="\t", names=cols)

except pd.errors.EmptyDataError:

    df = pd.DataFrame(columns=cols)



# drop self-hits

df = df[df["q"] != df["s"]]



if df.empty:

    print("[INFO] No BLAST hits; all ORFs will be Z clusters.")

    parent = {i: i for i in all_ids}

else:

    # convert to numeric

    for col in ["pident", "alen", "qlen", "slen", "evalue"]:

        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["pident", "alen", "qlen", "slen", "evalue"])



    # coverage = alen / min(qlen, slen)

    minlen = df[["qlen", "slen"]].min(axis=1)

    df["cov"] = df["alen"] / minlen



    mask = (

        (df["pident"] >= IDENT) &

        (df["cov"]    >= COV) &

        (df["evalue"] <= EVAL)

    )

    df_f = df[mask]



    print(f"[INFO] Total BLAST hits: {len(df)}")

    print(f"[INFO] Passing filters: {len(df_f)}")



    # ---- union–find clustering ----

    parent = {}

    def find(x):

        parent.setdefault(x, x)

        if parent[x] != x:

            parent[x] = find(parent[x])

        return parent[x]



    def union(a, b):

        pa, pb = find(a), find(b)

        if pa != pb:

            parent[pb] = pa



    for _, r in df_f.iterrows():

        q = r["q"]

        s = r["s"]

        if q in meta and s in meta:

            union(q, s)



    # ensure singletons

    for ident in all_ids:

        parent.setdefault(ident, ident)



# ---- build clusters ----

clusters = defaultdict(list)

for ident in all_ids:

    # find root with path compression

    root = ident

    while parent[root] != root:

        root = parent[root]

    clusters[root].append(ident)



print(f"[INFO] Number of clusters: {len(clusters)}")



rows = []

R_count = 0

G_count = 0

Z_count = 0



for root, members in clusters.items():

    species_set = {meta[m][0] for m in members}

    # Count species and samples per species

    species_list = sorted(species_set)

    species_count = len(species_list)

    samples_per_species = defaultdict(set)

    for m in members:

        sp, sa, _ = meta[m]

        samples_per_species[sp].add(sa)



    # Classification:

    # R = >= 2 species

    # G = 1 species & >= 2 samples

    # Z = 1 species & 1 sample

    if species_count >= 2:

        R_count += 1

        label = f"R{R_count}"

    else:

        sp = species_list[0]

        n_samples = len(samples_per_species[sp])

        if n_samples >= 2:

            G_count += 1

            label = f"G{G_count}"

        else:

            Z_count += 1

            label = f"Z{Z_count}"



    for m in members:

        sp, sa, orf = meta[m]

        rows.append([sp, sa, orf, label])



out = pd.DataFrame(rows, columns=["species", "sample", "orf_id", "cluster_label"])

out.to_csv("unannotated_RGZ_clusters.tsv", sep="\t", index=False)



print("[DONE] Wrote: unannotated_RGZ_clusters.tsv")

print("R clusters:", R_count)

print("G clusters:", G_count)

print("Z clusters:", Z_count)

EOF



echo

echo "[ALL DONE] Completed RGZ clustering pipeline."


