#!/usr/bin/env python3

import os

import re

import glob

import csv

from pathlib import Path



# ========= CONFIG =========

BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/dn_ds/M8_vs_M8a"

OUT_PARAMS = os.path.join(BASE, "M8_gene_params.csv")

OUT_PROP   = os.path.join(BASE, "M8_prop_sites.csv")



# Matches:

#   kappa (ts/tv) =  2.39469

RE_KAPPA = re.compile(r"kappa\s*\(ts/tv\)\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)")



#   Parameters in M8 (beta&w>1):

#     p0 =   0.99999  p =   0.35702 q =  42.92683

#    (p1 =   0.00001) w =   1.00000

RE_P0PQ = re.compile(

    r"p0\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+"

    r"p\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+"

    r"q\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"

)

RE_P1W = re.compile(

    r"\(p1\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\)\s*"

    r"w\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"

)



#   MLEs of dN/dS (w) for site classes (K=11)

#   p: 0.1 0.1 ... ; w: 0.00000 ...

RE_P_LINE = re.compile(r"^\s*p:\s*(.+)$")

RE_W_LINE = re.compile(r"^\s*w:\s*(.+)$")



def parse_float_list(s: str):

    # split on whitespace, convert to float, ignore empty pieces

    vals = []

    for tok in s.strip().split():

        # Strip possible trailing commas or other oddities

        tok = tok.strip(",")

        if tok:

            vals.append(float(tok))

    return vals



def parse_mlc(path):

    """

    Returns:

      {

        'kappa': float or None,

        'p0': float or None,

        'p': float or None,

        'q': float or None,

        'p1': float or None,

        'w': float or None,

        'p_vec': [floats] or None,

        'w_vec': [floats] or None,

      }

    """

    with open(path, "r", errors="ignore") as f:

        text = f.read()



    # Extract kappa

    kappa = None

    m = RE_KAPPA.search(text)

    if m:

        kappa = float(m.group(1))



    # Extract p0, p, q

    p0 = p = q = None

    m = RE_P0PQ.search(text)

    if m:

        p0, p, q = map(float, m.groups())



    # Extract p1 and w (the >1 class omega)

    p1 = w = None

    m = RE_P1W.search(text)

    if m:

        p1, w = map(float, m.groups())



    # Find p: and w: lines nearest to "MLEs of dN/dS..."

    p_vec = None

    w_vec = None

    # We’ll scan line by line to keep order and avoid false positives

    lines = text.splitlines()

    got_mles = False

    for i, line in enumerate(lines):

        if "MLEs of dN/dS (w) for site classes" in line:

            got_mles = True

            # next few lines should contain p: and w:

            # We’ll search forward a small window

            for j in range(i+1, min(i+10, len(lines))):

                lp = RE_P_LINE.search(lines[j])

                if lp:

                    p_vec = parse_float_list(lp.group(1))

                lw = RE_W_LINE.search(lines[j])

                if lw:

                    w_vec = parse_float_list(lw.group(1))

            break



    return {

        "kappa": kappa, "p0": p0, "p": p, "q": q, "p1": p1, "w": w,

        "p_vec": p_vec, "w_vec": w_vec

    }



def guess_gene_name(mlc_path: str) -> str:

    """

    Try to infer the gene name from the filename first (e.g., atp9_M8.mlc),

    else fall back to the parent folder name (e.g., .../atp9/M8/).

    """

    stem = Path(mlc_path).stem  # e.g., "atp9_M8"

    # If it ends with _M8, strip that

    if stem.endswith("_M8"):

        gene = stem[:-3]

    else:

        gene = stem

    if not gene:

        # fallback to parent of 'M8' dir

        parts = Path(mlc_path).parts

        try:

            idx = parts.index("M8")

            gene = parts[idx-1]

        except ValueError:

            gene = Path(mlc_path).parent.name

    return gene



def main():

    mlc_paths = glob.glob(os.path.join(BASE, "*", "M8", "*_M8.mlc"))

    if not mlc_paths:

        raise SystemExit(f"No *_M8.mlc files found under {BASE}/*/M8/")



    # Collect rows

    gene_param_rows = []

    prop_rows = []  # long format: gene, class, p, w



    for mlc in sorted(mlc_paths):

        gene = guess_gene_name(mlc)

        parsed = parse_mlc(mlc)



        gene_param_rows.append({

            "gene": gene,

            "kappa": parsed["kappa"],

            "p0": parsed["p0"],

            "p": parsed["p"],

            "q": parsed["q"],

            "p1": parsed["p1"],

            "w": parsed["w"]

        })



        p_vec = parsed["p_vec"] or []

        w_vec = parsed["w_vec"] or []

        K = max(len(p_vec), len(w_vec))

        # pad with Nones if lengths differ

        if len(p_vec) < K: p_vec = p_vec + [None]*(K-len(p_vec))

        if len(w_vec) < K: w_vec = w_vec + [None]*(K-len(w_vec))



        for k in range(K):

            prop_rows.append({

                "gene": gene,

                "site_class": k,

                "p": p_vec[k],

                "w": w_vec[k]

            })



    # Write outputs

    with open(OUT_PARAMS, "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=["gene","kappa","p0","p","q","p1","w"])

        writer.writeheader()

        writer.writerows(gene_param_rows)



    with open(OUT_PROP, "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=["gene","site_class","p","w"])

        writer.writeheader()

        writer.writerows(prop_rows)



    print(f"Wrote:\n  {OUT_PARAMS}\n  {OUT_PROP}")



if __name__ == "__main__":

    main()


