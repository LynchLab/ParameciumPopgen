import os

from Bio import SeqIO



# === INPUTS ===

BASE_DIR = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/protein_orf"

OUT_DIR  = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/phylogeny_with_outgroup/gene_fasta"

os.makedirs(OUT_DIR, exist_ok=True)



# Optional: lock to exactly these 46 genes (and this order)

GENE_LIST = None  # e.g., ["cox1", "cox2", "cob", "atp9", "nadh1-a", ...]



# Samples to exclude

EXCLUDE_SAMPLES = {

    "SC22-33__primaurelia",

    "WCL6-6__Plynn",

    "WL_5-8__triaurelia",

    "SWL_1-8__Pmulti",

    "LP7-7__triaurelia",

    "LP5A_2__triaurelia",

    "DP3B_1__biaurelia",

    "3_III_5__Plynn",
"M03___pfoki","SWL_1-3___multi","SWL_1-26___multi","SWL_1-29___caudatum"
}



def gene_from_filename(fname: str) -> str:

    """Files are '<sample>_<gene>.faa' (sample may contain underscores) → split from the right to get gene."""

    base = os.path.splitext(fname)[0]

    parts = base.rsplit("_", 1)

    return parts[-1] if len(parts) >= 2 else base



# 1) Discover samples and collect per-sample gene seqs

samples = []

sample_gene_seqs = {}  # sample -> {gene: seq_str}



for sample in sorted(os.listdir(BASE_DIR)):

    spath = os.path.join(BASE_DIR, sample)

    if not os.path.isdir(spath):

        continue

    if sample in EXCLUDE_SAMPLES:

        print(f"[SKIP SAMPLE] {sample}")

        continue



    samples.append(sample)

    gmap = {}

    for fname in os.listdir(spath):

        if not fname.endswith(".faa"):

            continue

        gene = gene_from_filename(fname)

        fpath = os.path.join(spath, fname)

        try:

            rec = next(SeqIO.parse(fpath, "fasta"))

        except StopIteration:

            continue

        # Uppercase, mask internal stops as X

        seq = str(rec.seq).upper().replace("*", "X")

        gmap[gene] = seq

    sample_gene_seqs[sample] = gmap



# 2) Determine gene set

if GENE_LIST is None:

    genes = sorted({g for gmap in sample_gene_seqs.values() for g in gmap.keys()})

else:

    genes = list(GENE_LIST)



print(f"Found {len(samples)} samples (after exclusion) and {len(genes)} genes to export.")



# 3) Write one FASTA per gene with all available samples

written = 0

for gene in genes:

    out_path = os.path.join(OUT_DIR, f"{gene}.faa")

    lines = []

    present = 0

    for sample in samples:

        seq = sample_gene_seqs.get(sample, {}).get(gene)

        if not seq:

            continue

        header = sample.replace(" ", "_")

        lines.append(f">{header}\n{seq}\n")

        present += 1



    if present > 0:

        with open(out_path, "w") as fh:

            fh.writelines(lines)

        written += 1

        print(f"[OK] {gene}: {present} sequences -> {out_path}")

    else:

        print(f"[SKIP] {gene}: no sequences found in any sample")



print(f"Done. Wrote {written} per-gene FASTA files to: {OUT_DIR}")


