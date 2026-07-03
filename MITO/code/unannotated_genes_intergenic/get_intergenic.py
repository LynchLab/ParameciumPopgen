#!/usr/bin/env python3

import pandas as pd

from pathlib import Path



ROOT = Path(

    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"

    "Ravina/additonal_species_identification/intergenic_region"

)



COORD_DIR = ROOT / "gene_coords_all"

OUTDIR = ROOT / "intergenic_tables"

OUTDIR.mkdir(exist_ok=True)





# -----------------------------

# CLEAN SAMPLE NAMES

# -----------------------------

def fix_sample_name(name):

    return (

        name.replace(".trnascan", "")

            .replace("_orf", "")

            .replace("__orf", "")

            .strip()

    )





# -----------------------------

# BUILD INTERGENIC TABLE

# -----------------------------

def get_intergenic_regions(sdf):

    """

    Given a dataframe of genes for a sample, return list of intergenic regions.

    """

    sdf = sdf.copy()

    sdf["start"] = sdf["start"].astype(int)

    sdf["end"] = sdf["end"].astype(int)

    sdf = sdf.sort_values("start").reset_index(drop=True)



    regions = []



    for g1, g2 in zip(sdf.iloc[:-1].itertuples(), sdf.iloc[1:].itertuples()):

        s = g1.end + 1

        e = g2.start - 1

        if s <= e:

            regions.append({

                "left_gene": g1.gene_name,

                "right_gene": g2.gene_name,

                "intergenic_name": f"{g1.gene_name}__{g2.gene_name}",

                "start": s,

                "end": e,

                "length": e - s + 1

            })



    return regions





# -----------------------------

# MAIN

# -----------------------------

def main():

    coord_files = sorted(COORD_DIR.glob("*_all_genes_annot_unannot_rrna_trna.csv"))



    if not coord_files:

        print("[ERROR] No gene coordinate files found.")

        return



    for csv in coord_files:

        species = csv.name.replace("_all_genes_annot_unannot_rrna_trna.csv", "")

        print(f"\n=== Processing species: {species} ===")



        df = pd.read_csv(csv, dtype=str)

        df["start"] = df["start"].astype(int)

        df["end"] = df["end"].astype(int)



        rows = []



        for sample_raw, sdf in df.groupby("sample"):

            sample_fixed = fix_sample_name(sample_raw)



            ig_regions = get_intergenic_regions(sdf)



            for r in ig_regions:

                rows.append({

                    "species": species,

                    "sample": sample_fixed,

                    **r

                })



        out_tsv = OUTDIR / f"{species}_intergenic_regions.tsv"

        pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)



        print(f"  → saved: {out_tsv}")



    print("\n✅ All intergenic TSV files created!")





if __name__ == "__main__":

    main()


