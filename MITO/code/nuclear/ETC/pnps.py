#!/usr/bin/env python3

"""

filter_pnps_by_ortholog_tsv_ids.py



Goal:

  From pcaud_paralogs_orthologs_transcript.tsv, collect ALL transcript IDs:

    - PCAU_protein column

    - every ID in orthologs_paralogs (comma-separated)

  Then, for each per-species pn/ps file in gene_pnps/gene_pnps/*.csv,

  keep only rows where gene is in that collected ID set.



Inputs:

  TSV:

    /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/pcaud_paralogs_orthologs_transcript.tsv



  pn/ps CSVs:

    /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/gene_pnps/gene_pnps/*.csv



Outputs:

  Writes filtered pn/ps CSVs into:

    /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/etc_pnps/

  Example:

    biaurelia.filtered_by_pcaud_orthologs.csv

    bursaria.filtered_by_pcaud_orthologs.csv

    ...

"""



import sys

from pathlib import Path

import pandas as pd





# -----------------------------

# CONFIG

# -----------------------------

TSV_IN = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/pcaud_paralogs_orthologs_transcript.tsv")



PNPS_DIR = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/gene_pnps/gene_pnps")



OUT_DIR = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/etc/etc_pnps")

OUT_DIR.mkdir(parents=True, exist_ok=True)



SUFFIX = "filtered_by_pcaud_orthologs"





def split_ids(cell) -> list[str]:

    """Split comma-separated IDs; safe on NaN."""

    if cell is None:

        return []

    s = str(cell).strip()

    if not s or s.lower() == "nan":

        return []

    return [x.strip() for x in s.split(",") if x.strip()]





def collect_ids_from_tsv(tsv_path: Path) -> set[str]:

    """Collect all transcript IDs from PCAU_protein and orthologs_paralogs fields."""

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)



    needed = {"PCAU_protein", "orthologs_paralogs"}

    missing = needed - set(df.columns)

    if missing:

        raise ValueError(f"TSV missing required columns: {sorted(missing)}")



    ids = set()



    # PCAU_protein column

    ids.update(x.strip() for x in df["PCAU_protein"].dropna().astype(str).tolist() if x.strip())



    # orthologs_paralogs list

    for cell in df["orthologs_paralogs"].tolist():

        ids.update(split_ids(cell))



    # remove junk like "nan"

    ids = {x for x in ids if x and x.lower() != "nan"}



    return ids





def main():

    if not TSV_IN.exists():

        sys.stderr.write(f"ERROR: TSV not found: {TSV_IN}\n")

        sys.exit(1)

    if not PNPS_DIR.exists():

        sys.stderr.write(f"ERROR: PNPS_DIR not found: {PNPS_DIR}\n")

        sys.exit(1)



    # 1) Collect IDs from ortholog TSV

    idset = collect_ids_from_tsv(TSV_IN)

    print(f"Collected {len(idset):,} unique IDs from ortholog TSV")



    # 2) Filter each pn/ps CSV

    pnps_files = sorted(PNPS_DIR.glob("*.csv"))

    if not pnps_files:

        sys.stderr.write(f"ERROR: No pn/ps CSV files found in: {PNPS_DIR}\n")

        sys.exit(1)



    for f in pnps_files:

        df = pd.read_csv(f, dtype=str)  # keep as strings to avoid float formatting issues

        if "gene" not in df.columns:

            sys.stderr.write(f"WARNING: Skipping {f.name} (no 'gene' column)\n")

            continue



        # keep original numeric columns as-is; filter by gene

        df["gene"] = df["gene"].astype(str).str.strip()

        keep = df["gene"].isin(idset)

        out = df.loc[keep].copy()



        out_path = OUT_DIR / f"{f.stem}.{SUFFIX}.csv"

        out.to_csv(out_path, index=False)



        print(f"{f.name}: kept {len(out):,} / {len(df):,}  ->  {out_path.name}")



    print(f"\n✅ All filtered files written to: {OUT_DIR}")





if __name__ == "__main__":

    main()


