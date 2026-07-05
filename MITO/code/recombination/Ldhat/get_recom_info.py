#!/usr/bin/env python3

import os

import re

import pandas as pd



base_dir = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination/gene_based"



# --- Regex Patterns ---

re_nseq = re.compile(r"Found\s+(\d+)\s+sequences")

re_length = re.compile(r"length\s+(\d+)")

re_div = re.compile(r"Estimated diversity.*?([\d\.]+)%")

re_ninfo = re.compile(r"Found\s+(\d+)\s+informative")

re_poly = re.compile(r"Window size is\s+(\d+)\s+polymorphic")

re_window = re.compile(r"Using a window size of\s+(\d+)")

re_k = re.compile(r"k as\s+(\d+)")

re_nss = re.compile(r"NSS:\s+([\d\.eE\+\-]+)")

re_maxchi = re.compile(r"Max Chi\^2:\s+([\d\.eE\+\-]+)")

re_phi_perm = re.compile(r"PHI \(Permutation\):\s+([\d\.eE\+\-]+)")

re_phi_norm = re.compile(r"PHI \(Normal\):\s+([\d\.eE\+\-]+)")



def extract(pattern, text):

    m = pattern.search(text)

    return m.group(1) if m else "NA"



rows = []



for species in os.listdir(base_dir):

    sp_dir = os.path.join(base_dir, species)

    if not os.path.isdir(sp_dir):

        continue



    for fname in os.listdir(sp_dir):

        if not fname.endswith("_phi_w100.log"):

            continue



        gene = fname.replace("_phi_w100.log", "")

        fpath = os.path.join(sp_dir, fname)



        with open(fpath, "r") as f:

            text = f.read()



        # ABSENT CASE

        if "fewer than 2 sequences" in text.lower():

            rows.append([

                species, gene,

                "NA","NA","NA","NA","NA","NA","NA",  # base stats

                "NA","NA","NA","NA"                  # recombination stats

            ])

            continue



        # Extract all stats

        window = extract(re_window, text)

        nseq = extract(re_nseq, text)

        nsites = extract(re_length, text)

        div = extract(re_div, text)

        ninfo = extract(re_ninfo, text)

        poly = extract(re_poly, text)

        k = extract(re_k, text)



        nss = extract(re_nss, text)

        maxchi = extract(re_maxchi, text)

        phi_perm = extract(re_phi_perm, text)

        phi_norm = extract(re_phi_norm, text)



        rows.append([

            species, gene,

            window, nseq, nsites, div, ninfo, poly, k,

            nss, maxchi, phi_perm, phi_norm

        ])



# Build DataFrame

df = pd.DataFrame(rows, columns=[

    "Species", "Gene",

    "Window", "nseq", "nsites", "div", "ninfo", "poly", "k",

    "NSS", "MaxChi2", "PHI_perm", "PHI_norm"

])



# Save CSV

outfile = os.path.join(base_dir, "phi_full_summary.csv")

df.to_csv(outfile, index=False)



print(f"✅ Saved summary CSV: {outfile}")


