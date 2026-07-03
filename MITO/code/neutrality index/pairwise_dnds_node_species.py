#!/usr/bin/env python3

"""

pairwise_dnds_node_species.py



Pairwise codon-change classification (dN/dS-style counts) between each species

and its closest internal node from closest_internal_node_per_species_mito_dnds.csv.



Behavior:

- filtered_sites/<gene>.txt is an EXCLUDE list per species:

    * Positions listed for that species are DROPPED from analysis (not counted in dN/dS).

    * These sites are still listed in the main TSV as Classification="Ignored",

      Ignored_Reason="ExcludedByList".

- Bad codons (gaps/N/?/non-ACGT or translating to stop/?) are also Ignored with reason "BadCodon".

- Out-of-range indices listed in the file are Ignored with reason "OutOfRange".

- If sequences aren’t a multiple of 3, trailing nts are trimmed and an Ignored row with

  Ignored_Reason="PostTrimTail" is added (no specific codon index).



Expected layout (run from ancestral_dnds_high_prob/):

  ./closest_internal_node_per_species_mito_dnds.csv

  ./filtered_sites/<gene>.txt                     # EXCLUDE list per species (1-based codon indices)

  ./<gene>/<gene>_all_sequences_cleaned_with_gaps.fasta      (preferred)

  (fallbacks also supported)



Outputs per gene (into <gene>/internal_node_dnds/):

  - <species>_vs_<node>_dnds.tsv            # includes analyzed rows + Ignored rows (with reasons)

  - <species>_vs_<node>_REMOVED_codons.tsv  # audit of ignored positions

Global summary:

  - summary_dnds_internal_nodes.tsv

"""



import os

import re

import pandas as pd

from Bio import SeqIO

from Bio.Seq import Seq

from Bio.Data import CodonTable



# === Config ===

CODON_TABLE_ID = 4  # Ciliate mitochondrial code (TGA=Trp; stops TAA/TAG)

MAPPING_FILE   = "closest_internal_node_per_species_mito_dnds.csv"

SUMMARY_OUT    = "summary_dnds_internal_nodes.tsv"

FILTER_DIR     = "filtered_sites"



codon_table = CodonTable.unambiguous_dna_by_id[CODON_TABLE_ID]

_species_key_re = re.compile(r"^[A-Za-z]+$")



# ---------- Helpers ----------



def load_mapping(path):

    """Return {gene: [(species, node_label), ...]} mapping."""

    df = pd.read_csv(path)

    df.columns = [c.strip().lower() for c in df.columns]

    if "internal_node" not in df.columns and "node" in df.columns:

        df.rename(columns={"node": "internal_node"}, inplace=True)

    # clean species like "5_prim" -> "prim"

    df["species"] = df["species"].astype(str).str.replace(r"^\d+_", "", regex=True)



    # normalize node labels to 'node_XX'

    def norm_node(x):

        x = str(x).strip()

        return x if x.startswith("node_") else f"node_{x}"

    df["internal_node"] = df["internal_node"].map(norm_node)



    mapping = {}

    for gene, grp in df.groupby("gene"):

        mapping[gene] = grp[["species", "internal_node"]].values.tolist()

    return mapping





def first_existing(*paths):

    for p in paths:

        if p and os.path.isfile(p):

            return p

    return None





def find_gene_fasta(gene):

    """

    Prefer <gene>/<gene>_all_sequences_cleaned_with_gaps.fasta

    Fallbacks:

      <gene>/all_sequences_cleaned_with_gaps.fasta

      <gene>/ancestral_MLbased/<gene>_all_sequences_cleaned_with_gaps.fasta

      <gene>/ancestral_MLbased/all_sequences_cleaned_with_gaps.fasta

    """

    preferred = os.path.join(gene, f"{gene}_all_sequences_cleaned_with_gaps.fasta")

    generic   = os.path.join(gene, "all_sequences_cleaned_with_gaps.fasta")

    ml_pref   = os.path.join(gene, "ancestral_MLbased", f"{gene}_all_sequences_cleaned_with_gaps.fasta")

    ml_gen    = os.path.join(gene, "ancestral_MLbased", "all_sequences_cleaned_with_gaps.fasta")

    return first_existing(preferred, generic, ml_pref, ml_gen)





def bad_codon(c):

    """True if unusable (gaps/N/?/non-ACGT or translates to stop/?)."""

    if len(c) != 3:

        return True

    u = c.upper()

    if any(x in u for x in "-N?"):

        return True

    if any(ch not in "ACGT" for ch in u):

        return True

    try:

        aa_ = str(Seq(u).translate(table=CODON_TABLE_ID))

    except Exception:

        return True

    return aa_ in {"*", "?"}





def aa(c):

    try:

        return str(Seq(c).translate(table=CODON_TABLE_ID))

    except Exception:

        return "?"





def load_species_filter_positions(gene):

    """

    Parse filtered_sites/<gene>.txt into {species: set(int codon_idx)}.



    File format (1-based codon indices separated by '-'; wrapped lines ok):

        prim 1-3-4-5-6-...

        pent 1-2-3-...

    Continuation lines (no species token) are appended to the current species.

    """

    path = os.path.join(FILTER_DIR, f"{gene}.txt")

    if not os.path.isfile(path):

        return {}  # No filter file => exclude none



    per = {}

    cur_species, acc = None, []



    with open(path) as fh:

        for raw in fh:

            line = raw.strip()

            if not line:

                continue

            parts = line.split()

            head = parts[0] if parts else ""

            # New species section?

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





# ---------- Main ----------



def main():

    if not os.path.isfile(MAPPING_FILE):

        raise SystemExit(f"Mapping file not found: {MAPPING_FILE}")



    gene_to_pairs = load_mapping(MAPPING_FILE)

    summary_rows = []



    # iterate gene subfolders in CWD, skipping non-gene dirs

    skip_dirs = {FILTER_DIR, "anc_site_probs"}

    genes = sorted(

        d for d in os.listdir(".")

        if os.path.isdir(d) and not d.startswith(".") and d not in skip_dirs

    )



    for gene in genes:

        if gene not in gene_to_pairs:

            continue



        fasta_path = find_gene_fasta(gene)

        if fasta_path is None:

            print(f"❌ Skipping {gene}: FASTA file not found.")

            continue



        sequences = {rec.id.strip(): str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")}

        species_filters = load_species_filter_positions(gene)



        out_dir = os.path.join(gene, "internal_node_dnds")

        os.makedirs(out_dir, exist_ok=True)



        for species, node in gene_to_pairs[gene]:

            seq1 = sequences.get(species)

            seq2 = sequences.get(node)

            if not seq1 or not seq2:

                print(f"❌ Missing sequence for {species} or {node} in {gene}")

                continue



            # --- Length bookkeeping ---

            prev_nt_len = min(len(seq1), len(seq2))

            trimmed_nt_len = prev_nt_len - (prev_nt_len % 3)

            s1 = seq1[:trimmed_nt_len]

            s2 = seq2[:trimmed_nt_len]



            total_codons = trimmed_nt_len // 3



            # EXCLUDE list for this species:

            # positions listed for the species are DROPPED (not included in dN/dS)

            exclude_set = species_filters.get(species, set())



            dn = ds = 0

            combined_rows = []  # analyzed + ignored rows for main TSV

            removed_rows  = []  # audit of ignored positions



            # If we trimmed tail nucleotides to achieve multiple-of-3, note it

            if prev_nt_len % 3 != 0:

                removed_rows.append(("n/a", "PostTrimTail"))



            # Walk codons (1-based index)

            for codon_idx in range(1, total_codons + 1):

                i = (codon_idx - 1) * 3

                c1 = s1[i:i+3]

                c2 = s2[i:i+3]

                pos_label = f"{i+1}-{i+3}"



                # Decide whether this site is ignored

                if codon_idx in exclude_set:

                    # Excluded by list

                    removed_rows.append((codon_idx, "ExcludedByList"))

                    a1 = aa(c1) if len(c1) == 3 and "?" not in c1 else "NA"

                    a2 = aa(c2) if len(c2) == 3 and "?" not in c2 else "NA"

                    combined_rows.append([pos_label, c1, c2, a1, a2, "Ignored", "ExcludedByList", codon_idx])

                    continue



                # Quality / biological filters

                if bad_codon(c1) or bad_codon(c2):

                    removed_rows.append((codon_idx, "BadCodon"))

                    a1 = aa(c1) if len(c1) == 3 and "?" not in c1 else "NA"

                    a2 = aa(c2) if len(c2) == 3 and "?" not in c2 else "NA"

                    combined_rows.append([pos_label, c1, c2, a1, a2, "Ignored", "BadCodon", codon_idx])

                    continue



                # Analyze

                a1 = aa(c1)

                a2 = aa(c2)

                if c1 == c2:

                    cls = "No change"

                elif a1 == a2:

                    cls = "Synonymous"; ds += 1

                else:

                    cls = "Non-synonymous"; dn += 1



                combined_rows.append([pos_label, c1, c2, a1, a2, cls, "", codon_idx])



            # Any exclude indices beyond available codons -> OutOfRange

            for pos in sorted(exclude_set):

                if pos < 1 or pos > total_codons:

                    removed_rows.append((pos, "OutOfRange"))

                    combined_rows.append(["NA", "NA", "NA", "NA", "NA", "Ignored", "OutOfRange", pos])



            # Optional explicit row for PostTrimTail in main file

            if any(r == "PostTrimTail" for _, r in removed_rows):

                combined_rows.append(["NA", "NA", "NA", "NA", "NA", "Ignored", "PostTrimTail", "n/a"])



            # Sort rows by codon index (numeric first, then 'n/a')

            def sort_key(row):

                idx = row[-1]

                if isinstance(idx, int):

                    return (0, idx)

                return (1, float('inf'))

            combined_rows.sort(key=sort_key)



            # --- Write per-comparison outputs ---

            base = f"{species}_vs_{node}"

            out_tsv = os.path.join(out_dir, f"{base}_dnds.tsv")

            pd.DataFrame(

                combined_rows,

                columns=[

                    "Codon_NT_Pos", f"Codon_{species}", f"Codon_{node}",

                    f"AA_{species}", f"AA_{node}",

                    "Classification", "Ignored_Reason", "Codon_Index"

                ]

            ).to_csv(out_tsv, sep="\t", index=False)



            removed_tsv = os.path.join(out_dir, f"{base}_REMOVED_codons.tsv")

            pd.DataFrame(

                removed_rows, columns=["Codon_Index", "Reason"]

            ).to_csv(removed_tsv, sep="\t", index=False)



            # --- Summary bookkeeping ---

            codons_included = sum(1 for row in combined_rows if row[5] in {"No change", "Synonymous", "Non-synonymous"})

            removed_excluded_count = sum(1 for _, r in removed_rows if r == "ExcludedByList")

            removed_bad_count      = sum(1 for _, r in removed_rows if r == "BadCodon")

            removed_oor_count      = sum(1 for _, r in removed_rows if r == "OutOfRange")

            removed_tail_flag      = int(any(r == "PostTrimTail" for _, r in removed_rows))

            nt_len_after_included  = codons_included * 3

            used_filter_file       = int(len(exclude_set) > 0)

            filter_file_path       = os.path.join(FILTER_DIR, f"{gene}.txt") if used_filter_file else ""



            summary_rows.append([

                gene, species, node,

                dn, ds,

                removed_excluded_count, removed_bad_count, removed_oor_count, removed_tail_flag,

                prev_nt_len, trimmed_nt_len,

                total_codons, codons_included, nt_len_after_included,

                used_filter_file, filter_file_path,

                out_tsv, removed_tsv

            ])



    # --- Global summary ---

    if summary_rows:

        pd.DataFrame(

            summary_rows,

            columns=[

                "Gene", "Species", "Node",

                "dN", "dS",

                "Sites_Ignored_ExcludedByList", "Sites_Ignored_Bad", "Sites_Ignored_OutOfRange", "Had_PostTrimTail",

                "Prev_NT_Len", "Trimmed_NT_Len",

                "Codons_Total", "Codons_Included", "NT_Len_After_Included",

                "Used_Filter_File", "Filter_File_Path",

                "Included_TSV", "Removed_TSV"

            ]

        ).to_csv(SUMMARY_OUT, sep="\t", index=False)

        print(f"✅ All genes processed. Summary saved as '{SUMMARY_OUT}'")

    else:

        print("⚠️ No comparisons written (check mapping/FASTA).")





if __name__ == "__main__":

    main()


