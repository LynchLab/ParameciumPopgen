#!/usr/bin/env python3

import pandas as pd

from pathlib import Path



# ======================

# CONFIG

# ======================



ROOT = Path(

    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"

    "Ravina/additonal_species_identification/intergenic_region"

)



INF_ROOT = Path(

    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"

    "Ravina/additonal_species_identification/infernal_rrna"

)



TRNA_ROOT = Path(

    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"

    "Ravina/additonal_species_identification/tRNA_genes"

)



UNANNOT_CLUSTER_FILE = ROOT / "RG_clustering" / "unannotated_RGZ_clusters.tsv"



OUTDIR = ROOT / "gene_coords_all"

OUTDIR.mkdir(exist_ok=True)





# ===============================================================

# Helper — resolve messy sample names in unannotated_RGZ_clusters

# ===============================================================



def resolve_sample_name(raw_sample, valid_samples):

    """

    Resolve a messy 'sample' string from unannotated_RGZ_clusters.tsv

    to an actual sample directory name in valid_samples.



    Strategy:

      1. If the cleaned string itself is in valid_samples, use it.

      2. If any valid sample is a substring of the cleaned string,

         use the longest matching one.

      3. Else, look for tokens that look like '*__*_orf' and see

         if any are in valid_samples; otherwise keep the longest such token.

      4. Else, fall back to first token (usually gives a warning later).

    """

    text = str(raw_sample).replace("\n", " ").strip().strip('"')



    # 1) direct match

    if text in valid_samples:

        return text



    # 2) any valid sample appears as substring

    substring_matches = [s for s in valid_samples if s in text]

    if substring_matches:

        return max(substring_matches, key=len)



    # 3) tokens that look like sample names (*__*_orf)

    tokens = text.split()

    candidates = [t.strip('"') for t in tokens if "__" in t and t.endswith("_orf")]



    # prefer any candidate that is a valid sample

    for c in candidates:

        if c in valid_samples:

            return c



    if candidates:

        return max(candidates, key=len)



    # 4) fallback: first token

    return text.split()[0] if text else text





# ===============================================================

# Annotated CDS from protein_annotated/*_blastp_with_coords.tsv

# ===============================================================



def parse_annotated(species_dir: Path, species: str):

    rows = []

    for sample_dir in sorted(species_dir.iterdir()):

        if not sample_dir.is_dir():

            continue

        sample_name = sample_dir.name



        pa_dir = sample_dir / "protein_annotated"

        if not pa_dir.is_dir():

            continue



        tsv_files = list(pa_dir.glob("*_blastp_with_coords.tsv"))

        if not tsv_files:

            continue

        ann = tsv_files[0]



        try:

            df = pd.read_csv(ann, sep="\t", dtype=str)

        except Exception as e:

            print(f"[WARN] Cannot read annotated TSV {ann}: {e}")

            continue



        required = {"query", "subject", "start", "stop", "align_len"}

        if not required.issubset(df.columns):

            print(f"[WARN] {ann} missing required columns {required - set(df.columns)}")

            continue



        for _, r in df.iterrows():

            subject = r["subject"]

            # e.g. PBIAU.V1-4.MITO.1.P040_cox1 -> 'cox1'

            gene_name = subject.split("_")[-1] if "_" in subject else subject



            try:

                start = int(r["start"])

                end = int(r["stop"])

                aa_len = int(r["align_len"])

            except Exception:

                print(f"[WARN] Non-integer coords/align_len in {ann} row: {r.to_dict()}")

                continue



            rows.append(

                {

                    "species": species,

                    "sample": sample_name,

                    "gene_name": gene_name,

                    "start": min(start, end),

                    "end": max(start, end),

                    "aa_length": aa_len,

                    "strand": "+" if start <= end else "-",

                    "feature_type": "annotated_cds",

                }

            )



    return rows





# ===============================================================

# Unannotated ORFs that are in R/G/Z clusters

# ===============================================================



def parse_unannotated(species: str, species_dir: Path):

    rows = []



    if not UNANNOT_CLUSTER_FILE.is_file():

        print(f"[WARN] No RGZ cluster file: {UNANNOT_CLUSTER_FILE}")

        return rows



    df = pd.read_csv(UNANNOT_CLUSTER_FILE, sep="\t", dtype=str)

    df = df[df["species"] == species].copy()

    if df.empty:

        return rows



    # set of real sample dirs for this species, e.g. "DP1G_2__biaurelia_orf"

    valid_samples = {

        d.name for d in species_dir.iterdir()

        if d.is_dir() and d.name.endswith("_orf")

    }



    for _, r in df.iterrows():

        raw_sample = r["sample"]

        sample = resolve_sample_name(raw_sample, valid_samples)



        if sample not in valid_samples:

            print(

                f"[WARN] Could not resolve valid sample for {species}: "

                f"raw='{raw_sample}' -> '{sample}'"

            )

            continue



        unann_f = species_dir / sample / "protein_unannotated" / "unannotated.csv"

        if not unann_f.is_file():

            print(f"[WARN] Missing unannotated.csv for {species}/{sample}")

            continue



        # robust read: try tab, then comma

        try:

            u = pd.read_csv(unann_f, sep="\t", dtype=str)

        except Exception:

            try:

                u = pd.read_csv(unann_f, sep=",", dtype=str)

            except Exception as e:

                print(f"[WARN] Cannot read {unann_f}: {e}")

                continue



        if "orf_id" not in u.columns:

            print(

                f"[WARN] {unann_f} has no 'orf_id' column; "

                f"columns={list(u.columns)}; skipping"

            )

            continue



        needed = {"orf_id", "start", "stop", "aa_len"}

        if not needed.issubset(u.columns):

            print(

                f"[WARN] {unann_f} missing cols {needed - set(u.columns)}; skipping"

            )

            continue



        u.set_index("orf_id", inplace=True)



        orf_id = r["orf_id"]

        cluster = r["cluster_label"]



        if orf_id not in u.index:

            print(

                f"[WARN] {species}/{sample}: ORF {orf_id} not found in {unann_f}"

            )

            continue



        urow = u.loc[orf_id]

        try:

            start = int(urow["start"])

            end = int(urow["stop"])

            aa_len = int(urow["aa_len"])

        except Exception:

            print(

                f"[WARN] Non-integer start/stop/aa_len in {unann_f} for {orf_id}"

            )

            continue



        rows.append(

            {

                "species": species,

                "sample": sample,

                "gene_name": cluster,  # R1 / G5 / Z3 etc.

                "start": min(start, end),

                "end": max(start, end),

                "aa_length": aa_len,

                "strand": "+" if start <= end else "-",

                "feature_type": "unannotated_ORF_cluster",

            }

        )



    return rows





# ===============================================================

# mt-rRNA from infernal_rrna

# ===============================================================



def guess_rrna_species_dir(species: str) -> Path | None:

    cand1 = INF_ROOT / f"{species}_fasta"

    cand2 = INF_ROOT / species

    if cand1.is_dir():

        return cand1

    if cand2.is_dir():

        return cand2

    return None





def parse_rrna(species: str):

    """

    From infernal_rrna/<species>_fasta/<sample>/<sample>.mt_rRNA.csv

    Columns: species,sample,contig,family,rfam_id,start,end,strand,length,score,evalue

    """

    rows = []

    spdir = guess_rrna_species_dir(species)

    if spdir is None:

        return rows



    for sample_dir in sorted(spdir.iterdir()):

        if not sample_dir.is_dir():

            continue

        sample_name = sample_dir.name

        mt_csv = sample_dir / f"{sample_name}.mt_rRNA.csv"

        if not mt_csv.is_file():

            continue



        try:

            df = pd.read_csv(mt_csv, dtype=str)

        except Exception as e:

            print(f"[WARN] Could not read {mt_csv}: {e}")

            continue



        required = {"sample", "family", "start", "end", "strand"}

        if not required.issubset(df.columns):

            print(f"[WARN] {mt_csv} missing {required - set(df.columns)}")

            continue



        for _, r in df.iterrows():

            try:

                start = int(r["start"])

                end = int(r["end"])

            except Exception:

                print(f"[WARN] Non-integer start/end in {mt_csv} row {r.to_dict()}")

                continue



            rows.append(

                {

                    "species": species,

                    "sample": r["sample"],

                    "gene_name": r["family"],  # e.g. LSU_rRNA_bacteria

                    "start": min(start, end),

                    "end": max(start, end),

                    "aa_length": "",

                    "strand": r["strand"],

                    "feature_type": "mt_rRNA",

                }

            )



    return rows





# ===============================================================

# tRNA from tRNAscan-SE (tRNA_genes)

# ===============================================================



def guess_trna_species_dir(species: str) -> Path | None:

    cand1 = TRNA_ROOT / f"{species}_fasta"

    cand2 = TRNA_ROOT / species

    if cand1.is_dir():

        return cand1

    if cand2.is_dir():

        return cand2

    return None





def parse_trna(species: str):

    """

    From tRNA_genes/<species>_fasta/<sample>.trnascan.tsv



    Example (after 3 header lines):

    Sequence  tRNA#  Begin  End   Type  Codon  IntronBegin  IntronEnd  Score  Note

    mito_biau_14 1   5743   5814  Phe   GAA   0            0          81.5   -

    """

    rows = []

    spdir = guess_trna_species_dir(species)

    if spdir is None:

        return rows



    for tsv in sorted(spdir.glob("*.trnascan.tsv")):

        sample = tsv.stem  # e.g. biaurelia_V1_4



        try:

            with open(tsv) as fh:

                # skip 3 header lines

                for _ in range(3):

                    next(fh, None)



                for line in fh:

                    line = line.strip()

                    if not line:

                        continue

                    cols = line.split()

                    if len(cols) < 8:

                        continue



                    try:

                        start = int(cols[2])

                        end = int(cols[3])

                    except ValueError:

                        continue



                    trna_type = cols[4]  # e.g. Phe, Sup, Tyr

                    gene_name = f"tRNA-{trna_type}"



                    rows.append(

                        {

                            "species": species,

                            "sample": sample,

                            "gene_name": gene_name,

                            "start": min(start, end),

                            "end": max(start, end),

                            "aa_length": "",

                            "strand": "+",  # strand not in TSV; assume +

                            "feature_type": "tRNA",

                        }

                    )

        except Exception as e:

            print(f"[WARN] Could not parse {tsv}: {e}")



    return rows





# ===============================================================

# MAIN

# ===============================================================



def main():

    exclude = {

        "RG_clustering",

        "unannotated_analysis",

        "gene_coords_all",

    }

    species_dirs = [

        p

        for p in sorted(ROOT.iterdir())

        if p.is_dir() and p.name not in exclude and not p.name.startswith(".")

    ]



    for spdir in species_dirs:

        species = spdir.name

        print(f"\n=== Processing {species} ===")



        rows = []

        rows.extend(parse_annotated(spdir, species))

        rows.extend(parse_unannotated(species, spdir))

        rows.extend(parse_rrna(species))

        rows.extend(parse_trna(species))



        if not rows:

            print(f"[WARN] No data collected for {species}")

            continue



        df = pd.DataFrame(rows)

        df = df.sort_values(["sample", "start", "end"], ascending=[True, True, True])



        outfile = OUTDIR / f"{species}_all_genes_annot_unannot_rrna_trna.csv"

        df.to_csv(outfile, index=False)

        print(f"  → Saved: {outfile}")



    print("\nAll species processed.")





if __name__ == "__main__":

    main()


