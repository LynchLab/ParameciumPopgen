#!/usr/bin/env python3

import os
import csv
from Bio import SeqIO

# -------------------------
# Paths
# -------------------------
CSV_DIR   = "prot_id"
FASTA_DIR = "../../oma_input"
OUT_DIR   = "protein_seq"

os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------
# CSV → FASTA mapping
# -------------------------
csv_to_fasta = {
    "biaurelia.csv": "biau.fa",
    "quindecaurelia.csv": "biau.fa",

    "multimicronucleatum.csv": "mult.fa",
    "fokini.csv": "mult.fa",

    "bursaria.csv": "burs.fa",
    "tribursaria.csv": "burs.fa",

    "primaurelia.csv": "prim_Ir4.fa",
    "triaurelia.csv": "prim_Ir4.fa",

    "caudatum.csv": "caud.fa",
    "sexaurelia.csv": "sexa.fa",
    "tetraurelia.csv": "tetr.fa",
}

# -------------------------
# Load FASTA files once
# -------------------------
fasta_cache = {}

def load_fasta(fa_path):
    if fa_path not in fasta_cache:
        fasta_cache[fa_path] = SeqIO.to_dict(SeqIO.parse(fa_path, "fasta"))
    return fasta_cache[fa_path]

def safe_fname(s: str) -> str:
    # gene IDs usually are safe, but just in case:
    return "".join(c if c.isalnum() or c in "._-+" else "_" for c in s)

# -------------------------
# Process each CSV
# -------------------------
for csv_file, fasta_file in csv_to_fasta.items():
    species = os.path.splitext(csv_file)[0]  # e.g., biaurelia
    species_dir = os.path.join(OUT_DIR, species)
    os.makedirs(species_dir, exist_ok=True)

    csv_path   = os.path.join(CSV_DIR, csv_file)
    fasta_path = os.path.join(FASTA_DIR, fasta_file)

    if not os.path.exists(csv_path):
        print(f"⚠️  Missing CSV: {csv_path}")
        continue

    if not os.path.exists(fasta_path):
        print(f"⚠️  Missing FASTA: {fasta_path}")
        continue

    records = load_fasta(fasta_path)

    found = 0
    missing = 0

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if "gene" not in reader.fieldnames:
            print(f"❌ {csv_file}: missing required header 'gene' (found {reader.fieldnames})")
            continue

        for row in reader:
            gene_id = row["gene"].strip()
            if not gene_id:
                continue

            if gene_id not in records:
                missing += 1
                continue

            rec = records[gene_id]
            out_path = os.path.join(species_dir, f"{safe_fname(gene_id)}.fa")
            with open(out_path, "w") as out:
                SeqIO.write(rec, out, "fasta")
            found += 1

    print(f"✅ {species}: wrote {found} protein FASTAs to {species_dir} (missing: {missing})")

