#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO

# ============================================================
# FULL PATHS
# ============================================================

INPUT_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/"
    "dnds_analysis/input/4_codon_alignment"
)

OUT_CSV = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/"
    "dnds_analysis/input/oma_species.csv"
)

# ============================================================
# ONLY ALLOWED TAXA (ORDER MATTERS)
# ============================================================

ALLOWED = ["PBIA", "PBUR", "PCAU", "PMMN", "PPRIM", "PSEX", "PTET"]

# ============================================================
# MAIN
# ============================================================

lines = ["oma_id,species"]

for fasta in sorted(INPUT_DIR.glob("OMA*.codon.aln.fna")):
    oma_id = fasta.stem.replace(".codon.aln", "")

    present = set()
    for rec in SeqIO.parse(fasta, "fasta"):
        if rec.id in ALLOWED:
            present.add(rec.id)

    if not present:
        continue

    ordered = [x for x in ALLOWED if x in present]
    lines.append(",".join([oma_id] + ordered))

OUT_CSV.write_text("\n".join(lines) + "\n")

print(f"✅ Wrote {len(lines)-1} rows to {OUT_CSV}")

