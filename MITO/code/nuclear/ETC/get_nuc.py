#!/usr/bin/env python3

import os
import csv
from Bio import SeqIO

# -------------------------
# Paths (run from dnds/etc)
# -------------------------
CSV_DIR   = "id"
FASTA_DIR = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/reference_assemblies/CDS"
OUT_DIR   = "nucleotide_seq"

os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------
# CSV -> CDS FASTA mapping (your rules)
# -------------------------
csv_to_fasta = {
    "biaurelia.csv": "biau.fa",
    "quindecaurelia.csv": "biau.fa",

    "caudatum.csv": "caud.fa",

    "multimicronucleatum.csv": "mult.fa",
    "fokini.csv": "mult.fa",

    "bursaria.csv": "burs.fa",
    "tribursaria.csv": "burs.fa",

    "tetraurelia.csv": "tetr.fa",

    "primaurelia.csv": "prim.fa",
    "triaurelia.csv": "prim.fa",

    "sexaurelia.csv": "sexa.fa",
}

# -------------------------
# FASTA cache
# -------------------------
fasta_cache = {}

def load_fasta(fa_path):
    """Load fasta into dict keyed by record.id once per file."""
    if fa_path not in fasta_cache:
        fasta_cache[fa_path] = SeqIO.to_dict(SeqIO.parse(fa_path, "fasta"))
    return fasta_cache[fa_path]

def safe_fname(s: str) -> str:
    """Make a filesystem-safe filename."""
    return "".join(c if c.isalnum() or c in "._-+" else "_" for c in s)

def read_gene_ids(csv_path: str):
    """Read 'gene' column, with a fallback for 1-column/no-header files."""
    ids = []
    with open(csv_path, newline="") as f:
        # Try DictReader first
        reader = csv.DictReader(f)
        if reader.fieldnames and "gene" in reader.fieldnames:
            for row in reader:
                gid = (row.get("gene") or "").strip()
                if gid:
                    ids.append(gid)
            return ids

    # Fallback: treat as one ID per line (possibly with header)
    with open(csv_path) as f2:
        for line in f2:
            s = line.strip()
            if not s or s.lower() == "gene":
                continue
            ids.append(s)
    return ids

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
    gene_ids = read_gene_ids(csv_path)

    found = 0
    missing = 0

    for gene_id in gene_ids:
        if gene_id not in records:
            missing += 1
            continue

        rec = records[gene_id]
        out_path = os.path.join(species_dir, f"{safe_fname(gene_id)}.fa")
        with open(out_path, "w") as out:
            SeqIO.write(rec, out, "fasta")
        found += 1

    print(f"✅ {species}: wrote {found} nucleotide FASTAs → {species_dir} (missing: {missing})")

