#!/usr/bin/env python3
"""
gc_content_fasta.py
Compute GC content (%) for one or more FASTA files.

Usage:
  python gc_content_fasta.py input.fasta
  python gc_content_fasta.py *.fa *.fasta
"""

from pathlib import Path
import sys

VALID = set("ACGTacgt")

def gc_percent(seq: str) -> float:
    seq = "".join([b for b in seq if b in VALID])  # keep only A/C/G/T
    if not seq:
        return float("nan")
    seq_u = seq.upper()
    gc = seq_u.count("G") + seq_u.count("C")
    return 100.0 * gc / len(seq_u)

def read_fasta(path: Path):
    header = None
    chunks = []
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)

def main():
    if len(sys.argv) < 2:
        print("Usage: python gc_content_fasta.py <fasta1> [fasta2 ...]")
        sys.exit(1)

    paths = [Path(p) for p in sys.argv[1:]]
    for p in paths:
        if not p.exists():
            print(f"ERROR: not found: {p}", file=sys.stderr)
            continue

        total_len = 0
        total_gc = 0

        print(f"\n# File: {p}")
        print("seq_id\tlength\tgc_percent")

        for seq_id, seq in read_fasta(p):
            seq_clean = [b for b in seq if b in VALID]
            L = len(seq_clean)
            if L == 0:
                gc = float("nan")
                print(f"{seq_id}\t0\tNaN")
                continue

            seq_u = "".join(seq_clean).upper()
            gc_count = seq_u.count("G") + seq_u.count("C")
            gc = 100.0 * gc_count / L

            total_len += L
            total_gc += gc_count

            print(f"{seq_id}\t{L}\t{gc:.3f}")

        if total_len > 0:
            file_gc = 100.0 * total_gc / total_len
            print(f"# Overall\t{total_len}\t{file_gc:.3f}")
        else:
            print("# Overall\t0\tNaN")

if __name__ == "__main__":
    main()
