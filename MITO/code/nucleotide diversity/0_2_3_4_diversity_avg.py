#!/usr/bin/env python3

import os

import sys

import pandas as pd





def summarize_one_species(species_dir: str):

    """

    Given a species directory like pn_ps_new/biaurelia,

    read per_site_pi/<species>_per_site_pi_ALL_GENES.csv

    and compute site-weighted pi for 1-,2-,3-,4-fold sites.

    """

    species = os.path.basename(os.path.abspath(species_dir))

    per_site_dir = os.path.join(species_dir, "per_site_pi")

    csv_path = os.path.join(per_site_dir, f"{species}_per_site_pi_ALL_GENES.csv")



    if not os.path.exists(csv_path):

        print(f"⚠️  Skipping {species}: missing {csv_path}")

        return None



    print(f"🔍 Processing {species}: {csv_path}")

    df = pd.read_csv(csv_path)



    # Drop sites where pi_site couldn't be computed

    df = df.dropna(subset=["pi_site"])



    out = {"species": species}



    # For each degeneracy class, take mean of pi_site over all sites of that class

    for d, col in [(1, "is_1fold"),

                   (2, "is_2fold"),

                   (3, "is_3fold"),

                   (4, "is_4fold")]:

        mask = df[col] == 1

        n_sites = int(mask.sum())



        if n_sites > 0:

            mean_pi = df.loc[mask, "pi_site"].mean()

        else:

            mean_pi = float("nan")



        out[f"n_{d}fold"] = n_sites

        out[f"pi_{d}fold"] = mean_pi



    return out





def main():

    # Root dir with species folders (biaurelia, caudatum, ...)

    root_dir = sys.argv[1] if len(sys.argv) > 1 else "."



    rows = []



    # Loop over everything in pn_ps_new/, pick things that look like species dirs

    for name in sorted(os.listdir(root_dir)):

        species_dir = os.path.join(root_dir, name)

        if not os.path.isdir(species_dir):

            continue



        # Only treat it as a species if it has per_site_pi/ inside

        per_site_dir = os.path.join(species_dir, "per_site_pi")

        if not os.path.isdir(per_site_dir):

            continue



        row = summarize_one_species(species_dir)

        if row is not None:

            rows.append(row)



    if not rows:

        print("❌ No species processed. Check directory structure.")

        sys.exit(1)



    out_df = pd.DataFrame(rows)

    out_csv = os.path.join(root_dir, "ALL_SPECIES_weighted_pi_by_degeneracy.csv")

    out_df.to_csv(out_csv, index=False)



    print(f"✅ Saved species-level weighted π table: {out_csv}")





if __name__ == "__main__":

    main()


