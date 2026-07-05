#!/usr/bin/env python3

"""

get_syn_nonsyn_sites_3.py



Compute synonymous / non-synonymous site *opportunities* for each

species ↔ closest internal node pair, using the SAME included sites as

pairwise_dnds_node_species.py in EXCLUDE mode.



Key behavior to match pairwise script:

- Use mitochondrial genetic code 4 (Ciliates).

- Trim to min(species_len, node_len), then enforce multiple-of-3.

- EXCLUDE list from filtered_sites/<gene>.txt applied by the species key

  (e.g., "18_pent" -> key "pent"). Positions in the list are DROPPED.

- A site is ignored if EITHER codon (species or node) is "bad":

  contains -, N, ?, non-ACGT, or translates to stop/? under code 4.

- Count site opportunities FROM THE NODE codon only, at included sites.



Outputs per gene (in <gene>/anc_dn_ds_cal_3/):

  - sites_<species>_vs_<node>.csv         (per-position opportunities, included sites only)

  - ignored_<species>_vs_<node>.tsv       (ignored codon indices + reason)

  - summary_sites_by_pair.tsv             (one row per species↔node pair)



Run from: ancestral_dnds_high_prob/

"""



import os

import re

import pandas as pd

from Bio import SeqIO

from Bio.Seq import Seq

from Bio.Data import CodonTable



# =========================

# Config

# =========================

CODON_TABLE_ID = 4  # Ciliate mitochondrial code (NCBI Table 4)

MAPPING_FILE   = "closest_internal_node_per_species_mito_dnds.csv"

FILTER_DIR     = "filtered_sites"



# =========================

# Helpers

# =========================

_species_key_re = re.compile(r"^[A-Za-z]+$")



def first_existing(*paths):

    for p in paths:

        if p and os.path.isfile(p):

            return p

    return None



def find_gene_fasta(gene):

    preferred = os.path.join(gene, f"{gene}_all_sequences_cleaned_with_gaps.fasta")

    generic   = os.path.join(gene, "all_sequences_cleaned_with_gaps.fasta")

    ml_pref   = os.path.join(gene, "ancestral_MLbased", f"{gene}_all_sequences_cleaned_with_gaps.fasta")

    ml_gen    = os.path.join(gene, "ancestral_MLbased", "all_sequences_cleaned_with_gaps.fasta")

    return first_existing(preferred, generic, ml_pref, ml_gen)



def load_mapping(path):

    """Return {gene: [(species_label, species_key, node_label), ...]}.



    species_label is used to try exact FASTA match first.

    species_key = stripped leading digits+underscore (e.g., '18_pent' -> 'pent') for filters and fallback lookup.

    node_label is normalized to 'node_XX'.

    """

    df = pd.read_csv(path)

    df.columns = [c.strip().lower() for c in df.columns]

    if "internal_node" not in df.columns and "node" in df.columns:

        df = df.rename(columns={"node": "internal_node"})



    def canon_species_key(s):

        s = str(s).strip()

        return re.sub(r"^\d+_", "", s)  # "18_pent" -> "pent"



    df["species_key"] = df["species"].map(canon_species_key)

    df["internal_node"] = df["internal_node"].astype(str).apply(

        lambda x: x if x.startswith("node_") else f"node_{x}"

    )



    out = {}

    for gene, grp in df.groupby("gene"):

        out[gene] = [(row["species"], row["species_key"], row["internal_node"]) for _, row in grp.iterrows()]

    return out



def load_species_filter_positions(gene):

    """filtered_sites/<gene>.txt  -> {species_key: set(int codon_idx)}  (EXCLUDE list)."""

    path = os.path.join(FILTER_DIR, f"{gene}.txt")

    if not os.path.isfile(path):

        return {}

    per = {}

    cur_species, acc = None, []

    with open(path) as fh:

        for raw in fh:

            line = raw.strip()

            if not line:

                continue

            parts = line.split()

            head = parts[0] if parts else ""

            if head and _species_key_re.match(head) and len(parts) >= 2:

                if cur_species and acc:

                    joined = "-".join(acc)

                    per[cur_species] = {

                        int(x) for x in joined.replace(" ", "").split("-") if x.isdigit()

                    }

                cur_species = head.strip()

                acc = [" ".join(parts[1:])]

            else:

                if cur_species:

                    acc.append(line)

    if cur_species and acc:

        joined = "-".join(acc)

        per[cur_species] = {

            int(x) for x in joined.replace(" ", "").split("-") if x.isdigit()

        }

    return per



# codon / AA helpers

NUCS = "ACGT"

ct = CodonTable.unambiguous_dna_by_id[CODON_TABLE_ID]

forward = ct.forward_table

stops   = set(ct.stop_codons)



def translate_codon(c):

    """Return AA; '*' for stop; '?' for invalid."""

    u = c.upper()

    if len(u) != 3 or any(x in u for x in "-N?") or any(ch not in "ACGT" for ch in u):

        return "?"

    if u in stops:

        return "*"

    return forward.get(u, "?")



def bad_codon(c):

    aa = translate_codon(c)

    return aa in {"?", "*"}



def site_opportunities_from_codon(codon):

    """Return (syn_frac_total, nonsyn_frac_total, detail_rows_for_this_codon)."""

    aa0 = translate_codon(codon)

    syn_total = 0.0

    nonsyn_total = 0.0

    rows = []

    for pos in range(3):

        syn_count = 0

        muts, aas, changed = [], [], []

        for n in NUCS:

            if n == codon[pos]:

                continue

            alt = codon[:pos] + n + codon[pos+1:]

            aa = translate_codon(alt)

            muts.append(alt)

            aas.append(aa)

            if aa not in {"?", "*"} and aa != aa0:

                changed.append(aa)

            if aa == aa0:

                syn_count += 1

        nonsyn_count = 3 - syn_count

        syn_total += syn_count / 3.0

        nonsyn_total += nonsyn_count / 3.0

        rows.append({

            "Position in Codon": pos + 1,

            "Original Codon": codon,

            "Original AA": aa0,

            "Mutated Codons": "; ".join(muts),

            "Mutated AAs": ", ".join(aas),

            "Changed AAs": ", ".join(sorted(set(changed))),

            "Synonymous Fraction": f"{syn_count}/3 = {syn_count/3:.2f}",

            "Non-synonymous Fraction": f"{nonsyn_count}/3 = {nonsyn_count/3:.2f}",

        })

    return syn_total, nonsyn_total, rows



def build_seq_indices(records):

    """Return (seq_dict, lower_index) for flexible label resolution."""

    seqs = {}

    lower_index = {}

    for rec in records:

        sid = rec.id.strip()

        seqs[sid] = str(rec.seq).upper()

        lower_index[sid.lower()] = sid

    return seqs, lower_index



def resolve_seq(label, lower_index, seqs):

    """Try exact, case-insensitive, and digit-stripped species key."""

    cand = [label, label.lower(), re.sub(r"^\d+_", "", label), re.sub(r"^\d+_", "", label).lower()]

    for c in cand:

        if c in seqs:

            return seqs[c], c

        if c in lower_index:

            return seqs[lower_index[c]], lower_index[c]

    return None, None



# =========================

# Main

# =========================

def main():

    if not os.path.isfile(MAPPING_FILE):

        raise SystemExit(f"Mapping file not found: {MAPPING_FILE}")

    gene_to_pairs = load_mapping(MAPPING_FILE)



    # iterate gene subfolders (skip non-genes)

    skip = {FILTER_DIR, "anc_site_probs"}

    genes = sorted(d for d in os.listdir(".") if os.path.isdir(d) and d not in skip and not d.startswith("."))



    for gene in genes:

        if gene not in gene_to_pairs:

            continue



        fasta = find_gene_fasta(gene)

        if fasta is None:

            print(f"❌ {gene}: FASTA not found; skipping.")

            continue



        records = list(SeqIO.parse(fasta, "fasta"))

        seqs, lower_index = build_seq_indices(records)

        filters = load_species_filter_positions(gene)



        out_dir = os.path.join(gene, "anc_dn_ds_cal_3")

        os.makedirs(out_dir, exist_ok=True)

        per_pair_summary = []



        # For quick diagnostics

        avail_nodes = sorted([sid for sid in seqs if sid.lower().startswith("node_")])



        for species_label, species_key, node_label in gene_to_pairs[gene]:

            sp_seq, sp_id = resolve_seq(species_label, lower_index, seqs)

            nd_seq, nd_id = resolve_seq(node_label,   lower_index, seqs)



            if sp_seq is None or nd_seq is None:

                # Helpful message with hints

                missing = []

                if sp_seq is None:

                    missing.append(f"{species_label} (key: {species_key})")

                if nd_seq is None:

                    missing.append(f"{node_label}")

                hint = ""

                if nd_seq is None and avail_nodes:

                    hint = f" | nodes in FASTA: {', '.join(avail_nodes[:10])}" + (" ..." if len(avail_nodes) > 10 else "")

                print(f"⚠️ {gene}: missing sequence for " + " or ".join(missing) + f"; skipping pair.{hint}")

                continue



            # Trim to min length and enforce multiple-of-3 (like pairwise)

            prev_nt_len = min(len(sp_seq), len(nd_seq))

            trimmed_nt_len = prev_nt_len - (prev_nt_len % 3)

            sp_aln = sp_seq[:trimmed_nt_len]

            nd_aln = nd_seq[:trimmed_nt_len]

            total_codons = trimmed_nt_len // 3



            # EXCLUDE list by species_key

            exclude = filters.get(species_key, set())

            allowed_positions = set(range(1, total_codons + 1)) - set(exclude)



            # bookkeeping

            ignored_rows = []  # [Codon_Index, Reason]

            if prev_nt_len % 3 != 0:

                ignored_rows.append(["n/a", "PostTrimTail"])



            # Per-position opportunity rows (included sites only)

            detail_rows = []

            dS_total = 0.0

            dN_total = 0.0

            included_codons = 0



            # Walk codons in alignment coordinates

            for idx in range(1, total_codons + 1):

                i = (idx - 1) * 3

                c_sp = sp_aln[i:i+3]

                c_nd = nd_aln[i:i+3]



                # Exclusion first (match pairwise behavior)

                if idx not in allowed_positions:

                    ignored_rows.append([idx, "ExcludedByList"])

                    continue



                # Bad codon if EITHER sequence bad at this position (match pairwise)

                if bad_codon(c_sp) or bad_codon(c_nd):

                    ignored_rows.append([idx, "BadCodon"])

                    continue



                # Count opportunities from the NODE codon only

                syn_add, nonsyn_add, pos_rows = site_opportunities_from_codon(c_nd)

                for pr in pos_rows:

                    pr_out = {

                        "Gene": gene,

                        "Species": species_label,

                        "Species_ResolvedID": sp_id,

                        "Node": node_label,

                        "Node_ResolvedID": nd_id,

                        "Codon #": idx,

                        "Codon_NT_Pos": f"{i+1}-{i+3}",

                        **pr

                    }

                    detail_rows.append(pr_out)



                dS_total += syn_add

                dN_total += nonsyn_add

                included_codons += 1



            # Out-of-range requested excludes (if any)

            for pos in sorted(exclude):

                if pos < 1 or pos > total_codons:

                    ignored_rows.append([pos, "OutOfRange"])



            # Write per-pair files

            base = f"{species_label}_vs_{node_label}"

            safe_base = base.replace("/", "_")



            sites_csv = os.path.join(out_dir, f"sites_{safe_base}.csv")

            pd.DataFrame(detail_rows).to_csv(sites_csv, index=False)



            ignored_tsv = os.path.join(out_dir, f"ignored_{safe_base}.tsv")

            pd.DataFrame(ignored_rows, columns=["Codon_Index", "Reason"]).to_csv(ignored_tsv, sep="\t", index=False)



            ratio = (dN_total / dS_total) if dS_total > 0 else float("inf")

            per_pair_summary.append({

                "Gene": gene,

                "Species": species_label,

                "Species_ResolvedID": sp_id,

                "Species_FilterKey": species_key,

                "Node": node_label,

                "Node_ResolvedID": nd_id,

                "Synonymous Sites": round(dS_total, 3),

                "Non-synonymous Sites": round(dN_total, 3),

                "dN/dS (opportunity)": ("Inf" if ratio == float("inf") else round(ratio, 4)),

                "Prev_NT_Len": prev_nt_len,

                "Trimmed_NT_Len": trimmed_nt_len,

                "Codons_Total": total_codons,

                "Codons_Included": included_codons,

                "Ignored_ExcludedByList": sum(1 for _, r in ignored_rows if r == "ExcludedByList"),

                "Ignored_BadCodon":       sum(1 for _, r in ignored_rows if r == "BadCodon"),

                "Ignored_OutOfRange":     sum(1 for _, r in ignored_rows if r == "OutOfRange"),

                "Had_PostTrimTail":       int(any(r == "PostTrimTail" for _, r in ignored_rows)),

                "Sites_CSV": sites_csv,

                "Ignored_TSV": ignored_tsv,

            })



            print(f"✅ {gene}: {species_label} (→ {sp_id}) vs {node_label} (→ {nd_id}) -> {sites_csv}")



        # Per-gene summary

        if per_pair_summary:

            summary_csv = os.path.join(out_dir, "summary_sites_by_pair.tsv")

            pd.DataFrame(per_pair_summary).to_csv(summary_csv, sep="\t", index=False)

            print(f"📄 {gene}: summary written to {summary_csv}")

        else:

            print(f"⚠️ {gene}: no pairs summarized (check mapping/FASTA).")



if __name__ == "__main__":

    main()


