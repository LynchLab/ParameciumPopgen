#!/usr/bin/env python3

import os

import sys

import math

import csv

from pathlib import Path

from typing import List, Tuple



import pandas as pd

from Bio import AlignIO



# ----------------------------

# CONFIG: set this to the species folder that contains:

#   codon_alignment/  (FASTA files like cox2_codon_aligned.fna)

#   degeneracy_result/ (TXT files like cox2.txt with digits 1/2/3/4/0)

# ----------------------------

# Example:

# species_dir = "/data/.../pn_ps_new/biaurelia"

species_dir = sys.argv[1] if len(sys.argv) > 1 else "."



CODON_ALIGN_DIR = "codon_alignment"

DEGENERACY_DIR  = "degeneracy_result"

OUT_DIR         = "per_site_pi"

VALID_NUCS      = set("ATGC")  # ignore gaps/ambiguity (N, -, etc.)



os.makedirs(os.path.join(species_dir, OUT_DIR), exist_ok=True)



def per_site_pi_from_alignment(aln) -> List[float]:

    """

    Compute per-site nucleotide diversity (pi) for a multiple sequence alignment.

    pi(site) = (# of pairwise differences) / (# of valid pairs), considering only A/T/G/C.

    Gaps or ambiguous bases (N, -) are skipped.

    Returns a list of length = alignment length.

    """

    nseq = len(aln)

    ncol = aln.get_alignment_length()

    pis = []



    # Convert to list of strings for speed

    seqs = ["".join(rec) if isinstance(rec, list) else str(rec.seq) for rec in aln]



    for col in range(ncol):

        col_bases = [seqs[i][col].upper() for i in range(nseq)]

        # keep only valid nucleotides

        vals = [b for b in col_bases if b in VALID_NUCS]

        m = len(vals)

        if m < 2:

            pis.append(float("nan"))  # not enough data

            continue



        # Count pairwise differences

        diffs = 0

        total = 0

        for i in range(m):

            for j in range(i + 1, m):

                total += 1

                if vals[i] != vals[j]:

                    diffs += 1

        pi = diffs / total if total > 0 else float("nan")

        pis.append(pi)



    return pis



def read_degeneracy_vector(path: str) -> List[int]:

    """

    Read degeneracy string (e.g., '214214...') from TXT and return list of ints.

    We keep only characters in {0,1,2,3,4}. (0 can mean non-coding/invalid columns)

    """

    with open(path, "r") as f:

        s = f.read().strip()

    # Keep only digits

    digits = [c for c in s if c.isdigit()]

    # Map to ints

    return [int(c) for c in digits]



def one_hot_from_deg(d: int) -> Tuple[int, int, int, int]:

    """

    Given degeneracy class (1,2,3,4; anything else -> 0s), return one-hot (is1,is2,is3,is4).

    """

    return (

        1 if d == 1 else 0,

        1 if d == 2 else 0,

        1 if d == 3 else 0,

        1 if d == 4 else 0,

    )



def gene_name_from_fasta(fname: str) -> str:

    """

    Extract gene name from filename like 'cox2_codon_aligned.fna' -> 'cox2'

    """

    base = os.path.basename(fname)

    if base.endswith("_codon_aligned.fna"):

        return base.replace("_codon_aligned.fna", "")

    return os.path.splitext(base)[0]



def main():

    codon_dir = os.path.join(species_dir, CODON_ALIGN_DIR)

    degen_dir = os.path.join(species_dir, DEGENERACY_DIR)

    out_dir   = os.path.join(species_dir, OUT_DIR)



    if not os.path.isdir(codon_dir):

        sys.exit(f"❌ Missing directory: {codon_dir}")

    if not os.path.isdir(degen_dir):

        sys.exit(f"❌ Missing directory: {degen_dir}")



    fasta_files = sorted([

        os.path.join(codon_dir, f)

        for f in os.listdir(codon_dir)

        if f.endswith("_codon_aligned.fna")

    ])

    if not fasta_files:

        sys.exit(f"❌ No *_codon_aligned.fna files found in: {codon_dir}")



    all_rows = []

    species_name = os.path.basename(os.path.abspath(species_dir))



    for fpath in fasta_files:

        gene = gene_name_from_fasta(fpath)

        degen_path = os.path.join(degen_dir, f"{gene}.txt")

        if not os.path.exists(degen_path):

            print(f"⚠️  Skipping {gene}: no degeneracy file {degen_path}")

            continue



        try:

            aln = AlignIO.read(fpath, "fasta")

        except Exception as e:

            print(f"⚠️  Skipping {gene}: failed to read alignment ({e})")

            continue



        # Per-site pi

        pi_sites = per_site_pi_from_alignment(aln)

        aln_len = len(pi_sites)



        # Degeneracy vector

        deg_vec = read_degeneracy_vector(degen_path)



        # If lengths differ, align by the min length, but warn

        L = min(aln_len, len(deg_vec))

        if aln_len != len(deg_vec):

            print(f"⚠️  {gene}: alignment length ({aln_len}) != degeneracy length ({len(deg_vec)}). Truncating to {L}.")



        # Build per-site rows

        rows = []

        for i in range(L):

            d = deg_vec[i]

            is1, is2, is3, is4 = one_hot_from_deg(d)

            rows.append({

                "species": species_name,

                "gene": gene,

                "site": i + 1,  # 1-based

                "is_1fold": is1,

                "is_2fold": is2,

                "is_3fold": is3,

                "is_4fold": is4,

                "pi_site": pi_sites[i]

            })



        # Save per-gene CSV

        gene_df = pd.DataFrame(rows)

        out_csv = os.path.join(out_dir, f"{gene}_per_site_pi.csv")

        gene_df.to_csv(out_csv, index=False)

        print(f"✅ Saved per-site table: {out_csv}")



        # Accumulate for combined file

        all_rows.extend(rows)



    # Combined CSV for the species

    if all_rows:

        big_df = pd.DataFrame(all_rows)

        combined_csv = os.path.join(out_dir, f"{species_name}_per_site_pi_ALL_GENES.csv")

        big_df.to_csv(combined_csv, index=False)

        print(f"✅ Combined table saved: {combined_csv}")

    else:

        print("⚠️ No per-site rows generated (no overlaps between alignments and degeneracy files?).")



if __name__ == "__main__":

    main()


