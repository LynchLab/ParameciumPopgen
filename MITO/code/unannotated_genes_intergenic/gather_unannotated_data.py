#!/usr/bin/env python3

import os

import sys

import pandas as pd



def main():

    # base directory: where you run the script (intergenic_region)

    base_dir = os.path.abspath(".")



    # loop over species directories (e.g. biaurelia, primaurelia, ...)

    for species in sorted(os.listdir(base_dir)):

        species_dir = os.path.join(base_dir, species)

        if not os.path.isdir(species_dir):

            continue



        all_dfs = []



        # loop over sample directories inside each species

        for sample in sorted(os.listdir(species_dir)):

            sample_dir = os.path.join(species_dir, sample)

            if not os.path.isdir(sample_dir):

                continue



            unann_path = os.path.join(sample_dir, "protein_unannotated", "unannotated.csv")

            if not os.path.isfile(unann_path):

                continue



            try:

                df = pd.read_csv(unann_path, sep="\t")

            except Exception as e:

                print(f"[WARN] Could not read {unann_path}: {e}", file=sys.stderr)

                continue



            # add species and sample columns

            df.insert(0, "sample", sample)

            df.insert(0, "species", species)



            all_dfs.append(df)



        if not all_dfs:

            print(f"[INFO] No unannotated ORFs found for species {species}, skipping.")

            continue



        out_df = pd.concat(all_dfs, ignore_index=True)



        out_path = os.path.join(base_dir, f"{species}_unannotated_orfs.tsv")

        out_df.to_csv(out_path, sep="\t", index=False)



        print(f"[DONE] Wrote {out_path} with {len(out_df)} rows.")



if __name__ == "__main__":

    main()


