#!/usr/bin/env python3

"""

species_avg_pnps_vs_ps.py



Compute per-species averages of:

  - mean pN/pS (from split Group A)

  - mean pS     (from split Group B)



Then test correlation ACROSS species.



Inputs:

  BASE_ROOT/<species>/codon_alignment/*_codon_aligned.fna



Outputs (all under OUT_ROOT):

  per_species_means.tsv

  species_correlation.tsv

"""



from pathlib import Path

from itertools import combinations

from collections import Counter

import math

import pandas as pd

from Bio import SeqIO

from Bio.Data import CodonTable



# ================= CONFIG =================

BASE_ROOT = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pn_ps_new")

OUT_ROOT  = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/grouping_pnps_ps")



MIN_SAMPLES = 4

MIN_CONS_FRAC = 0.60

GENETIC_CODE_ID = 4

SEED = 7

CORR_METHOD = "spearman"   # or "pearson"

# ==========================================



ctab = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]

AA_TABLE = ctab.forward_table

STOP_CODONS = set(ctab.stop_codons)

START_CODONS = [c for c in ctab.start_codons if c != "TTA"] + ["GTA"]



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

        s1, s2 = str(r1.seq), str(r2.seq)

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



def split(recs, seed):

    import random

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



        for fna in (sdir / "codon_alignment").glob("*_codon_aligned.fna"):

            recs = list(SeqIO.parse(fna, "fasta"))

            if len(recs) < MIN_SAMPLES:

                continue



            cons = consensus_sequence(recs, MIN_CONS_FRAC)

            S_sites, N_sites = site_totals(cons)

            if S_sites <= 0 or N_sites <= 0:

                continue



            A, B = split(recs, SEED + hash((sp, fna.name)) % 10000)

            if len(A) < 2 or len(B) < 2:

                continue



            synA, nonsynA, pA = pairwise_counts(A)

            synB, _, pB = pairwise_counts(B)

            if pA == 0 or pB == 0:

                continue



            pS_A = (synA / pA) / S_sites

            pN_A = (nonsynA / pA) / N_sites

            if pS_A > 0:

                pnps_vals.append(pN_A / pS_A)



            ps_vals.append((synB / pB) / S_sites)



        rows.append({

            "Species": sp,

            "Genes_used": len(pnps_vals),

            "Mean_pN_pS": pd.Series(pnps_vals).mean(),

            "Mean_pS": pd.Series(ps_vals).mean()

        })



    df = pd.DataFrame(rows).dropna()

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


