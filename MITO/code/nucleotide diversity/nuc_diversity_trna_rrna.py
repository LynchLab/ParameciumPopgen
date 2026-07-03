#!/usr/bin/env python3

"""

pi_trna_rrna_by_species.py



Reads *_pi_per_site.tsv files in each species' mafft/pi/ directory,

and computes site-weighted mean π for tRNA vs rRNA genes.



Output: pi_trna_rrna_by_species.csv



Columns:

    Species, tRNA_pi, rRNA_pi

"""



import os

import glob

from pathlib import Path

import pandas as pd



BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification"

INF  = f"{BASE}/infernal_rrna"



def classify_gene_type_from_pi_file_name(tsv_path: Path):

    """

    Classify gene type using the pi-per-site TSV filename.



    Examples:

        tRNA_Phe_pi_per_site.tsv  -> tRNA

        tRNA_Trp_pi_per_site.tsv  -> tRNA

        LSU_rRNA_pi_per_site.tsv  -> rRNA

        SSU_rRNA_pi_per_site.tsv  -> rRNA

    """

    name = tsv_path.name  # full filename

    stem = tsv_path.stem  # e.g. LSU_rRNA_pi_per_site



    if stem.startswith("tRNA_"):

        return "tRNA"

    # treat anything with 'rRNA' in the name as rRNA, as long as it is not tRNA

    if "rrna" in stem.lower():

        return "rRNA"



    # Default to rRNA if you know only tRNAs and rRNAs are present

    return "rRNA"



def main():

    rows = []



    for species_dir in sorted(glob.glob(f"{INF}/*")):

        if not os.path.isdir(species_dir):

            continue



        species = os.path.basename(species_dir)

        pi_dir = Path(species_dir) / "mafft" / "pi"

        if not pi_dir.is_dir():

            print(f"  ⚠️ No mafft/pi directory for {species}, skipping.")

            continue



        print(f"\n=== Species (π): {species} ===")



        # gather all *_pi_per_site.tsv

        tsv_files = sorted(glob.glob(str(pi_dir / "*_pi_per_site.tsv")))

        if not tsv_files:

            print(f"  ⚠️ No *_pi_per_site.tsv files for {species}, skipping.")

            continue



        # site-weighted accumulators

        trna_pi_sum = 0.0

        trna_site_count = 0

        rrna_pi_sum = 0.0

        rrna_site_count = 0



        for tsv in tsv_files:

            tsv_path = Path(tsv)

            gene_type = classify_gene_type_from_pi_file_name(tsv_path)



            # read Pi column

            df = pd.read_csv(tsv_path, sep="\t")

            if "Pi" not in df.columns:

                print(f"    ⚠️ File {tsv_path} missing 'Pi' column, skipping.")

                continue



            pi_values = df["Pi"].dropna()

            site_count = len(pi_values)

            if site_count == 0:

                continue



            pi_sum = pi_values.sum()



            if gene_type == "tRNA":

                trna_pi_sum += pi_sum

                trna_site_count += site_count

            elif gene_type == "rRNA":

                rrna_pi_sum += pi_sum

                rrna_site_count += site_count



        # final site-weighted means

        trna_pi_mean = trna_pi_sum / trna_site_count if trna_site_count > 0 else None

        rrna_pi_mean = rrna_pi_sum / rrna_site_count if rrna_site_count > 0 else None



        rows.append({

            "Species": species,

            "tRNA_pi": trna_pi_mean,

            "rRNA_pi": rrna_pi_mean

        })



    if not rows:

        print("⚠️ No π rows generated.")

        return



    df = pd.DataFrame(rows)

    out_path = Path(INF) / "pi_trna_rrna_by_species.csv"

    df.to_csv(out_path, index=False)

    print(f"\n📦 π summary written: {out_path}")



if __name__ == "__main__":

    main()


