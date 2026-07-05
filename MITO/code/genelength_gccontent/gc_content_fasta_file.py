#!/usr/bin/env python3
import argparse, pathlib, re, sys
from collections import defaultdict

def read_fasta(path):
    seq = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq)

def gc_counts(seq):
    seq = seq.upper()
    # count only A/C/G/T; ignore Ns and anything else
    valid = re.findall(r'[ACGT]', seq)
    g = seq.count('G')
    c = seq.count('C')
    return len(valid), (g + c)

def main():
    p = argparse.ArgumentParser(description="Compute GC% for all genes and a specified subset.")
    p.add_argument("-d", "--dir", default=".", help="Directory containing *.fna (default: current dir)")
    p.add_argument("--subset", nargs="*", default=["rpL6", "nadh7", "ymf79", "rpS14"],
                   help="Gene basenames (case-insensitive) to treat as negative-strand subset")
    p.add_argument("-o", "--out_csv", default="per_gene_gc.csv", help="Per-gene CSV output")
    p.add_argument("-s", "--summary", default="gc_summary.txt", help="Summary text output")
    args = p.parse_args()

    d = pathlib.Path(args.dir)
    files = sorted(d.glob("*.fna"))
    if not files:
        print(f"[ERROR] No .fna files found in: {d}", file=sys.stderr)
        sys.exit(1)

    # Normalize subset names (just the gene token, case-insensitive)
    subset_set = {x.lower() for x in args.subset}

    per_gene = []  # (gene, length, gc_count, gc_pct, in_subset)
    total_len_all = 0
    total_gc_all = 0
    total_len_subset = 0
    total_gc_subset = 0

    for f in files:
        seq = read_fasta(f)
        length, gc = gc_counts(seq)

        # Infer gene name token from filename: take the last part before .fna
        # e.g., "11_II_7__sexaurelia_rpL6.fna" -> "rpL6"
        stem = f.stem
        # handle names that may contain underscores; take the trailing token after last "_"
        gene = stem.split("_")[-1]
        gene_lc = gene.lower()

        gc_pct = (gc / length * 100.0) if length > 0 else 0.0
        in_subset = (gene_lc in subset_set)

        per_gene.append((gene, length, gc, gc_pct, in_subset))

        total_len_all += length
        total_gc_all  += gc
        if in_subset:
            total_len_subset += length
            total_gc_subset  += gc

    all_gc_pct = (total_gc_all / total_len_all * 100.0) if total_len_all else 0.0
    sub_gc_pct = (total_gc_subset / total_len_subset * 100.0) if total_len_subset else 0.0

    # write per-gene CSV
    with open(args.out_csv, "w") as out:
        out.write("gene,length,gc_count,gc_pct,subset_negative\n")
        for gene, length, gc, gc_pct, in_subset in per_gene:
            out.write(f"{gene},{length},{gc},{gc_pct:.4f},{int(in_subset)}\n")

    # write summary
    with open(args.summary, "w") as out:
        out.write("GC Content Summary\n")
        out.write("==================\n")
        out.write(f"All genes combined:     length={total_len_all}, GC={total_gc_all}, GC%={all_gc_pct:.4f}\n")
        out.write(f"Negative-strand subset: length={total_len_subset}, GC={total_gc_subset}, GC%={sub_gc_pct:.4f}\n")
        out.write(f"Subset genes considered: {', '.join(sorted(args.subset))}\n")

    # also print concise summary to stdout
    print(f"[ALL]   length={total_len_all} GC={total_gc_all} GC%={all_gc_pct:.4f}")
    print(f"[NEG]   length={total_len_subset} GC={total_gc_subset} GC%={sub_gc_pct:.4f}")
    print(f"[INFO]  per-gene CSV -> {args.out_csv}")
    print(f"[INFO]  summary      -> {args.summary}")

if __name__ == "__main__":
    main()
