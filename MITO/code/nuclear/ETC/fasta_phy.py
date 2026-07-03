#!/usr/bin/env python3
"""
convert_oma_fasta_to_phylip.py

Convert all OMA*.codon.aln.fna FASTA alignments to sequential PHYLIP (.phy),
PAML-safe: strict equal-length check, simple names, one-line sequences.
"""

from __future__ import annotations
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = None
    seq_chunks: List[str] = []
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(seq_chunks).replace(" ", "").replace("\t", "").upper()
                    records.append((name, seq))
                name = line[1:].strip().split()[0]  # take first token only
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            seq = "".join(seq_chunks).replace(" ", "").replace("\t", "").upper()
            records.append((name, seq))
    return records


def sanitize_name(name: str, maxlen: int = 10) -> str:
    # Keep alnum/._- only; PHYLIP classic name field is 10 chars
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "", name)
    if not clean:
        clean = "seq"
    return clean[:maxlen]


def write_phylip_sequential(outpath: Path, records: List[Tuple[str, str]]) -> None:
    if len(records) < 2:
        raise ValueError("Fewer than 2 sequences; PHYLIP usually needs >=2.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        # show a helpful message with per-seq lengths
        details = ", ".join([f"{nm}:{len(seq)}" for nm, seq in records])
        raise ValueError(f"Length mismatch (must all be equal). {details}")

    L = next(iter(lengths))
    n = len(records)

    # Detect duplicated (sanitized) names which can confuse PAML
    names = [sanitize_name(nm) for nm, _ in records]
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate names after sanitization: {names}")

    with outpath.open("w") as out:
        out.write(f"{n} {L}\n")
        for (raw_name, seq), nm in zip(records, names):
            # PHYLIP sequential: name padded to 10 columns + space + sequence
            out.write(f"{nm.ljust(10)} {seq}\n")


def main() -> int:
    in_dir = Path(
        "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/"
        "nuclear_etc_ribo/dnds/etc_ribo/input_codon_alignment"
    )
    out_dir = Path(
        "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/"
        "nuclear_etc_ribo/dnds/etc_ribo/phylip_input"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(in_dir.glob("OMA*.codon.aln.fna"))
    if not fasta_files:
        print(f"❌ No files matching OMA*.codon.aln.fna in {in_dir}", file=sys.stderr)
        return 1

    n_ok = 0
    n_fail = 0

    for fa in fasta_files:
        oma = fa.name.replace(".codon.aln.fna", "")
        outpath = out_dir / f"{oma}.phy"
        try:
            recs = read_fasta(fa)
            write_phylip_sequential(outpath, recs)
            print(f"✅ {oma}: wrote {outpath}")
            n_ok += 1
        except Exception as e:
            print(f"❌ {oma}: FAILED ({e})", file=sys.stderr)
            n_fail += 1

    print(f"\nDone. Success: {n_ok}, Failed: {n_fail}. Output: {out_dir}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

