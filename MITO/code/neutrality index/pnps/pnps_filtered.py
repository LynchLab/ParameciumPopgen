#!/usr/bin/env python3

"""

pnps_filtered_3.py



Inputs (per species):

  <BASE_ROOT>/<species>/codon_alignment_filtered/*_codon_aligned_filtered.fna



Outputs (per species, written inside each species folder under BASE_ROOT):

  <BASE_ROOT>/<species>/pnps/consensus/<gene>_consensus.fna

  <BASE_ROOT>/<species>/pnps/sites_detail/<gene>_codon_mutation_table.csv

  <BASE_ROOT>/<species>/pnps/sites/<gene>_sites.csv

  <BASE_ROOT>/<species>/pnps/pairwise/<gene>_pairwise_summary.csv

  <BASE_ROOT>/<species>/pnps/pairwise_detail/<gene>/<sampleA_vs_sampleB>_codons.csv   (if SAVE_PAIRWISE_CODON_TABLES=True)

  <BASE_ROOT>/<species>/pnps/pairwise_across_genes/<species>_sample_pairwise_summary.csv

  <BASE_ROOT>/<species>/pnps/summaries/<species>_pnps_mean_summary.csv



Combined (all species):

  <BASE_ROOT>/pnps_mean_summary_all_species.csv

  <BASE_ROOT>/pnps_mean_summary_all_species_no_missing_ratio.csv

"""



from pathlib import Path

from itertools import combinations

from collections import Counter, defaultdict

import pandas as pd

from Bio import SeqIO

from Bio.Data import CodonTable



# ================= CONFIG =================

BASE_ROOT = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pnps_sites_0.6")

ALIGN_DIR_NAME = "codon_alignment_filtered"  # << use filtered alignments

ALIGN_GLOB = "*_codon_aligned_filtered.fna"  # << file name pattern

GENE_NAME_CUT = "_codon_aligned_filtered.fna"  # << strip to get gene name



MIN_CONS_FRAC = 0.6

GENETIC_CODE_ID = 4  # Mold/Protozoan (Paramecium mt)

SAVE_PAIRWISE_CODON_TABLES = False

# ==========================================



# --- genetic code + custom start/stop rules ---

ctab = CodonTable.unambiguous_dna_by_id[GENETIC_CODE_ID]

AA_TABLE = ctab.forward_table

STOP_CODONS = set(ctab.stop_codons)

# Custom start rules: remove TTA, add GTA

START_CODONS = [c for c in ctab.start_codons if c != "TTA"] + ["GTA"]





def chunks_of3(s: str):

    s = s.upper()

    return [s[i:i+3] for i in range(0, len(s), 3) if len(s[i:i+3]) == 3]





def aa_call(codon: str, is_start: bool):

    if len(codon) != 3 or "N" in codon or "-" in codon:

        return None

    if codon in STOP_CODONS:

        return "Stop"

    if is_start and codon in START_CODONS:

        return "M"

    return AA_TABLE.get(codon, None)





def write_fasta(path: Path, name: str, seq: str):

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as fh:

        fh.write(f">{name}\n{seq}\n")





def consensus_sequence(records, minc=0.6):

    if not records:

        return ""

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





def build_codon_site_table(cons_seq: str, species: str, gene: str):

    """Per-codon table + totals; start rules apply only to first KEPT codon."""

    nucs = "ACGT"

    S_total = 0.0

    N_total = 0.0

    rows = []

    kept_idx = 0



    codons_aln = chunks_of3(cons_seq.replace("-", ""))  # ignore alignment gaps for consensus evaluation

    for aln_idx, codon in enumerate(codons_aln, start=1):

        if "N" in codon or codon in STOP_CODONS or codon not in AA_TABLE:

            continue

        kept_idx += 1

        is_start = (kept_idx == 1)

        orig_aa = aa_call(codon, is_start=is_start)

        if orig_aa in (None, "Stop"):

            kept_idx -= 1

            continue



        for pos in range(3):

            mutated_codons, mutated_aas, changed_aas = [], [], []

            syn_count = 0

            for n in nucs:

                if n == codon[pos]:

                    continue

                alt = codon[:pos] + n + codon[pos+1:]

                aa = aa_call(alt, is_start=is_start)

                mutated_codons.append(alt)

                mutated_aas.append(aa if aa is not None else "NA")

                if aa == orig_aa:

                    syn_count += 1

                elif aa not in (None, "Stop"):

                    changed_aas.append(aa)

            nonsyn_count = 3 - syn_count

            S_total += syn_count / 3.0

            N_total += nonsyn_count / 3.0

            rows.append({

                "Codon #": kept_idx,

                "Alignment Codon #": aln_idx,

                "Position in Codon": pos + 1,

                "Original Codon": codon,

                "Original AA": orig_aa,

                "Mutated Codons": "; ".join(mutated_codons),

                "Mutated AAs": ", ".join(map(str, mutated_aas)),

                "Changed AAs": ", ".join(sorted(set(changed_aas))),

                "Synonymous Fraction": f"{syn_count}/3 = {syn_count/3:.2f}",

                "Non-synonymous Fraction": f"{nonsyn_count}/3 = {nonsyn_count/3:.2f}",

            })



    site_row = {

        "Species": species,

        "Gene": gene,

        "Synonymous Sites": round(S_total, 6),

        "Non-synonymous Sites": round(N_total, 6),

    }

    return rows, site_row, S_total, N_total





def pairwise_classify(rec1, rec2):

    """Counts per pair; start applies to first NON-IGNORED codon in the pair."""

    s1, s2 = str(rec1.seq).upper(), str(rec2.seq).upper()

    L = min(len(s1), len(s2))

    codons1, codons2 = chunks_of3(s1[:L]), chunks_of3(s2[:L])



    rows = []

    syn = nonsyn = ident = ign = tot = 0

    start_assigned = False



    for idx, (c1, c2) in enumerate(zip(codons1, codons2), start=1):

        tot += 1

        ignore_reason = ""

        if "-" in c1 or "-" in c2 or "N" in c1 or "N" in c2:

            ign += 1

            classification = "Ignored"

            ignore_reason = "gap_or_N"

            aa1 = aa2 = None

        else:

            is_start = False

            if not start_assigned:

                is_start = True

                start_assigned = True

            aa1 = aa_call(c1, is_start=is_start)

            aa2 = aa_call(c2, is_start=is_start)

            if aa1 in (None, "Stop") or aa2 in (None, "Stop"):

                ign += 1

                classification = "Ignored"

                ignore_reason = "stop_or_invalid"

            else:

                if c1 == c2:

                    ident += 1

                    classification = "Identical"

                elif aa1 == aa2:

                    syn += 1

                    classification = "Synonymous"

                else:

                    nonsyn += 1

                    classification = "Non-synonymous"



        if SAVE_PAIRWISE_CODON_TABLES:

            rows.append({

                "Codon # (alignment)": idx,

                "Seq1 Codon": c1,

                "Seq2 Codon": c2,

                "Seq1 AA": ("Stop" if aa1 == "Stop" else (aa1 if aa1 else "")),

                "Seq2 AA": ("Stop" if aa2 == "Stop" else (aa2 if aa2 else "")),

                "Classification": classification,

                "IgnoreReason": ignore_reason

            })



    return (syn, nonsyn, ident, ign, tot), rows





def ensure_outdirs(species_dir: Path):

    """Create <species>/pnps/* inside BASE_ROOT."""

    root = species_dir / "pnps"

    (root / "consensus").mkdir(parents=True, exist_ok=True)

    (root / "sites_detail").mkdir(parents=True, exist_ok=True)

    (root / "sites").mkdir(parents=True, exist_ok=True)

    (root / "pairwise").mkdir(parents=True, exist_ok=True)

    (root / "pairwise_detail").mkdir(parents=True, exist_ok=True)

    (root / "pairwise_across_genes").mkdir(parents=True, exist_ok=True)

    (root / "summaries").mkdir(parents=True, exist_ok=True)

    return root





def main():

    # species = subfolders of BASE_ROOT that contain ALIGN_DIR_NAME

    species_dirs = [p for p in sorted(BASE_ROOT.iterdir())

                    if p.is_dir() and (p / ALIGN_DIR_NAME).is_dir()]

    if not species_dirs:

        print(f"[ERROR] No species with {ALIGN_DIR_NAME} under {BASE_ROOT}")

        return



    combined_rows = []



    for sdir in species_dirs:

        sp = sdir.name

        ca_dir = sdir / ALIGN_DIR_NAME

        out_root = ensure_outdirs(sdir)



        print(f"\n🔬 {sp}")

        gene_files = sorted(ca_dir.glob(ALIGN_GLOB))

        if not gene_files:

            print(f"  ⚠️ {sp}: no {ALIGN_GLOB}")

            continue



        pair_agg = defaultdict(lambda: {

            "Synonymous": 0,

            "Non-synonymous": 0,

            "Identical": 0,

            "Ignored": 0,

            "Total Codons": 0,

            "GenesUsed": 0

        })

        species_summary_rows = []



        for fna in gene_files:

            gene = fna.name.replace(GENE_NAME_CUT, "")

            recs = list(SeqIO.parse(str(fna), "fasta"))

            if not recs:

                print(f"  ⚠️ {sp}:{gene} empty alignment")

                continue



            # consensus + site totals

            cons_seq = consensus_sequence(recs, minc=MIN_CONS_FRAC)

            write_fasta(out_root / "consensus" / f"{gene}_consensus.fna",

                        f"{gene}_consensus", cons_seq)



            codon_rows, site_row, S_sites, N_sites = build_codon_site_table(cons_seq, sp, gene)

            pd.DataFrame(codon_rows).to_csv(out_root / "sites_detail" / f"{gene}_codon_mutation_table.csv", index=False)

            pd.DataFrame([site_row]).to_csv(out_root / "sites" / f"{gene}_sites.csv", index=False)



            # pairwise per gene

            if len(recs) < 2:

                pd.DataFrame(columns=[

                    "Species", "Gene", "Pair", "Synonymous", "Non-synonymous", "Identical",

                    "Ignored", "TotalCodons", "Synonymous_Sites", "NonSyn_Sites", "pS", "pN", "pN_pS"

                ]).to_csv(out_root / "pairwise" / f"{gene}_pairwise_summary.csv", index=False)



                species_summary_rows.append({

                    "Species": sp, "Gene": gene, "Pairs": 0,

                    "Synonymous_Sites": S_sites, "NonSyn_Sites": N_sites,

                    "Mean_Ds": float('nan'), "Mean_Dn": float('nan'),

                    "Mean_pS": float('nan'), "Mean_pN": float('nan'),

                    "pN_over_pS": float('nan')

                })

                continue



            pair_rows = []

            for r1, r2 in combinations(recs, 2):

                id1 = r1.id.split()[0].replace("__", "_")

                id2 = r2.id.split()[0].replace("__", "_")

                key = tuple(sorted([id1, id2]))



                counts, codon_detail = pairwise_classify(r1, r2)

                syn, nonsyn, ident, ign, tot = counts



                pS = (syn / S_sites) if S_sites > 0 else float('nan')

                pN = (nonsyn / N_sites) if N_sites > 0 else float('nan')

                ratio = (pN / pS) if (pS and pS == pS) else float('inf')



                if SAVE_PAIRWISE_CODON_TABLES:

                    d = out_root / "pairwise_detail" / gene

                    d.mkdir(parents=True, exist_ok=True)

                    pd.DataFrame(codon_detail).to_csv(d / f"{key[0]}_vs_{key[1]}_codons.csv", index=False)



                pair_rows.append({

                    "Species": sp, "Gene": gene, "Pair": f"{key[0]}_vs_{key[1]}",

                    "Synonymous": syn, "Non-synonymous": nonsyn,

                    "Identical": ident, "Ignored": ign, "TotalCodons": tot,

                    "Synonymous_Sites": S_sites, "NonSyn_Sites": N_sites,

                    "pS": pS, "pN": pN, "pN_pS": ratio

                })



                # accumulate across genes for this sample pair

                pair_agg[key]["Synonymous"] += syn

                pair_agg[key]["Non-synonymous"] += nonsyn

                pair_agg[key]["Identical"] += ident

                pair_agg[key]["Ignored"] += ign

                pair_agg[key]["Total Codons"] += tot

                pair_agg[key]["GenesUsed"] += 1



            df_pairs = pd.DataFrame(pair_rows)

            df_pairs.to_csv(out_root / "pairwise" / f"{gene}_pairwise_summary.csv", index=False)



            # mean counts -> per-site -> ratio

            mean_ds = df_pairs["Synonymous"].mean() if not df_pairs.empty else float('nan')

            mean_dn = df_pairs["Non-synonymous"].mean() if not df_pairs.empty else float('nan')

            mean_pS = (mean_ds / S_sites) if (S_sites > 0 and pd.notna(mean_ds)) else float('nan')

            mean_pN = (mean_dn / N_sites) if (N_sites > 0 and pd.notna(mean_dn)) else float('nan')

            mean_ratio = (mean_pN / mean_pS) if (pd.notna(mean_pS) and mean_pS != 0) else float('inf')



            species_summary_rows.append({

                "Species": sp, "Gene": gene, "Pairs": len(df_pairs),

                "Synonymous_Sites": S_sites, "NonSyn_Sites": N_sites,

                "Mean_Ds": mean_ds, "Mean_Dn": mean_dn,

                "Mean_pS": mean_pS, "Mean_pN": mean_pN, "pN_over_pS": mean_ratio

            })



        # per-species summaries

        if species_summary_rows:

            (sdir / "pnps" / "summaries").mkdir(parents=True, exist_ok=True)

            sp_sum_path = sdir / "pnps" / "summaries" / f"{sp}_pnps_mean_summary.csv"

            pd.DataFrame(species_summary_rows).sort_values(["Species", "Gene"]).to_csv(sp_sum_path, index=False)

            combined_rows.extend(species_summary_rows)



        # sample↔sample aggregated across genes

        if pair_agg:

            out_rows = []

            for (id1, id2), v in sorted(pair_agg.items()):

                out_rows.append({"Species": sp, "SamplePair": f"{id1}_vs_{id2}", **v})

            (sdir / "pnps" / "pairwise_across_genes").mkdir(parents=True, exist_ok=True)

            pd.DataFrame(out_rows).to_csv(

                sdir / "pnps" / "pairwise_across_genes" / f"{sp}_sample_pairwise_summary.csv",

                index=False

            )

            print("  ✅ wrote sample↔sample across-genes summary")



    # combined all-species at BASE_ROOT

    if combined_rows:

        all_path = BASE_ROOT / "pnps_mean_summary_all_species.csv"

        pd.DataFrame(combined_rows).to_csv(all_path, index=False)



        filtered_path = BASE_ROOT / "pnps_mean_summary_all_species_no_missing_ratio.csv"

        pd.DataFrame(combined_rows).dropna(subset=["pN_over_pS"]).to_csv(filtered_path, index=False)



        print(f"\n📦 Combined written: {all_path}")

        print(f"🧹 Filtered (drop missing pN_over_pS) written: {filtered_path}")



    print("\n✅ All done. Inputs were read from:", BASE_ROOT)





if __name__ == "__main__":

    main()
