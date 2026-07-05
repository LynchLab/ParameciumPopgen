#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
from collections import Counter
import math
import hashlib
import random
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable

# ================= CONFIG =================
BASE_ROOT = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pn_ps_new")
OUT_ROOT  = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/grouping_pnps_ps")
MIN_SAMPLES = 4
MIN_CONS_FRAC = 0.60
GENETIC_CODE_ID = 4
CORR_METHOD = "spearman"

N_SPLITS = 50         # number of repeated random splits per gene
BASE_SEED = 12345     # fixed for reproducibility
# ==========================================

ctab = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]
AA_TABLE = ctab.forward_table
STOP_CODONS = set(ctab.stop_codons)
START_CODONS = [c for c in ctab.start_codons if c != "TTA"] + ["GTA"]


def mean_se(vals):
    vals = [v for v in vals if v is not None and np.isfinite(v)]
    n = len(vals)
    if n == 0:
        return math.nan, math.nan, 0
    mean_val = float(np.mean(vals))
    if n == 1:
        return mean_val, math.nan, 1
    se_val = float(np.std(vals, ddof=1) / math.sqrt(n))
    return mean_val, se_val, n


def stable_hash_int(text):
    return int(hashlib.md5(text.encode()).hexdigest(), 16)


def chunks_of3(s):
    return [s[i:i+3] for i in range(0, len(s), 3) if len(s[i:i+3]) == 3]


def aa_call(codon, is_start):
    if len(codon) != 3 or "N" in codon or "-" in codon:
        return None
    if codon in STOP_CODONS:
        return "Stop"
    if is_start and codon in START_CODONS:
        return "M"
    return AA_TABLE.get(codon, None)


def consensus_sequence(records, minc=0.6):
    L = len(records[0].seq)
    out = []
    for i in range(L):
        col = [str(r.seq)[i].upper() for r in records]
        bases = [b for b in col if b in "ACGT"]
        if not bases:
            out.append("N")
            continue
        base, count = Counter(bases).most_common(1)[0]
        out.append(base if count / len(col) >= minc else "N")
    return "".join(out)


def site_totals(cons_seq):
    S = N = 0.0
    kept = 0

    for codon in chunks_of3(cons_seq.replace("-", "")):
        if "N" in codon or codon in STOP_CODONS or codon not in AA_TABLE:
            continue

        kept += 1
        is_start = (kept == 1)
        aa0 = aa_call(codon, is_start)

        if aa0 in (None, "Stop"):
            kept -= 1
            continue

        for pos in range(3):
            syn = 0
            for n in "ACGT":
                if n == codon[pos]:
                    continue
                alt = codon[:pos] + n + codon[pos+1:]
                aa = aa_call(alt, is_start)
                if aa == aa0:
                    syn += 1
            S += syn / 3
            N += (3 - syn) / 3

    return S, N


def pairwise_counts(recs):
    syn = nonsyn = pairs = 0

    for r1, r2 in combinations(recs, 2):
        s1, s2 = str(r1.seq).upper(), str(r2.seq).upper()
        start_used = False
        ps = pn = 0

        for c1, c2 in zip(chunks_of3(s1), chunks_of3(s2)):
            if "-" in c1 or "-" in c2 or "N" in c1 or "N" in c2:
                continue

            is_start = False
            if not start_used:
                is_start = True
                start_used = True

            aa1 = aa_call(c1, is_start)
            aa2 = aa_call(c2, is_start)

            if aa1 in (None, "Stop") or aa2 in (None, "Stop"):
                continue

            if c1 != c2:
                if aa1 == aa2:
                    ps += 1
                else:
                    pn += 1

        syn += ps
        nonsyn += pn
        pairs += 1

    return syn, nonsyn, pairs


def pairwise_ps_values(recs, S_sites):
    values = []
    if S_sites <= 0:
        return values

    for r1, r2 in combinations(recs, 2):
        s1, s2 = str(r1.seq).upper(), str(r2.seq).upper()
        start_used = False
        ps = 0

        for c1, c2 in zip(chunks_of3(s1), chunks_of3(s2)):
            if "-" in c1 or "-" in c2 or "N" in c1 or "N" in c2:
                continue

            is_start = False
            if not start_used:
                is_start = True
                start_used = True

            aa1 = aa_call(c1, is_start)
            aa2 = aa_call(c2, is_start)

            if aa1 in (None, "Stop") or aa2 in (None, "Stop"):
                continue

            if c1 != c2 and aa1 == aa2:
                ps += 1

        values.append(ps / S_sites)

    return values


def split_random(recs, seed):
    r = list(recs)
    random.Random(seed).shuffle(r)
    nA = len(r) // 2
    return r[:nA], r[nA:]


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []

    for sdir in sorted(BASE_ROOT.iterdir()):
        if not (sdir / "codon_alignment").is_dir():
            continue

        sp = sdir.name
        pnps_vals = []
        ps_vals = []
        ps_pairwise_se_vals = []

        for fna in sorted((sdir / "codon_alignment").glob("*_codon_aligned.fna")):
            recs = list(SeqIO.parse(fna, "fasta"))

            if len(recs) < MIN_SAMPLES:
                continue

            cons = consensus_sequence(recs, MIN_CONS_FRAC)
            S_sites, N_sites = site_totals(cons)

            if S_sites <= 0 or N_sites <= 0:
                continue

            gene_pnps_reps = []
            gene_ps_reps = []
            gene_pairwise_ps_se_reps = []

            gene_seed_base = BASE_SEED + stable_hash_int(f"{sp}::{fna.name}") % 1000000

            for i in range(N_SPLITS):
                split_seed = gene_seed_base + i
                A, B = split_random(recs, split_seed)

                if len(A) < 2 or len(B) < 2:
                    continue

                # Group A: pN/pS
                synA, nonsynA, pA = pairwise_counts(A)
                if pA > 0:
                    pS_A = (synA / pA) / S_sites
                    pN_A = (nonsynA / pA) / N_sites
                    if pS_A > 0 and np.isfinite(pS_A) and np.isfinite(pN_A):
                        ratio = pN_A / pS_A
                        if np.isfinite(ratio):
                            gene_pnps_reps.append(ratio)

                # Group B: pS
                synB, _, pB = pairwise_counts(B)
                if pB > 0:
                    pS_B = (synB / pB) / S_sites
                    if np.isfinite(pS_B):
                        gene_ps_reps.append(pS_B)

                    gene_pairwise_ps = pairwise_ps_values(B, S_sites)
                    _, gene_ps_se, n_pairs = mean_se(gene_pairwise_ps)
                    if n_pairs >= 2 and np.isfinite(gene_ps_se):
                        gene_pairwise_ps_se_reps.append(gene_ps_se)

            # average across repeated splits for this gene
            gene_mean_pnps, _, n_pnps_rep = mean_se(gene_pnps_reps)
            gene_mean_ps, _, n_ps_rep = mean_se(gene_ps_reps)
            gene_mean_pair_ps_se, _, n_pair_ps_rep = mean_se(gene_pairwise_ps_se_reps)

            if n_pnps_rep > 0 and np.isfinite(gene_mean_pnps):
                pnps_vals.append(gene_mean_pnps)

            if n_ps_rep > 0 and np.isfinite(gene_mean_ps):
                ps_vals.append(gene_mean_ps)

            if n_pair_ps_rep > 0 and np.isfinite(gene_mean_pair_ps_se):
                ps_pairwise_se_vals.append(gene_mean_pair_ps_se)

        mean_pnps, se_pnps, n_pnps = mean_se(pnps_vals)
        mean_ps, se_ps, n_ps = mean_se(ps_vals)
        mean_ps_pair_se, _, n_pairse = mean_se(ps_pairwise_se_vals)

        rows.append({
            "Species": sp,
            "Genes_used_pN_pS": n_pnps,
            "Mean_pN_pS": mean_pnps,
            "SE_pN_pS": se_pnps,
            "Genes_used_pS": n_ps,
            "Mean_pS": mean_ps,
            "SE_pS": se_ps,
            "Mean_gene_pairwise_SE_pS": mean_ps_pair_se,
            "Genes_with_pairwise_SE_pS": n_pairse
        })

    df = pd.DataFrame(rows).dropna(subset=["Mean_pN_pS", "Mean_pS"])
    df.to_csv(OUT_ROOT / "per_species_means.tsv", sep="\t", index=False)

    rho = df["Mean_pN_pS"].corr(df["Mean_pS"], method=CORR_METHOD)
    pd.DataFrame([{
        "Correlation_type": CORR_METHOD,
        "rho(pN/pS , pS)": rho,
        "N_species": len(df)
    }]).to_csv(OUT_ROOT / "species_correlation.tsv", sep="\t", index=False)

    print("✅ Outputs written to", OUT_ROOT)


if __name__ == "__main__":
    main()