#!/usr/bin/env python3

import csv

from pathlib import Path

from collections import defaultdict

from statistics import mean



ROOT = Path(".")  # run inside nucleotide_orf/

OUT_SAMPLE  = "gc_per_samples.csv"

OUT_SPECIES = "gc_per_species.csv"



VALID = set("ATGC")



def gc_percent(seq: str) -> float:

    """GC% over A/T/G/C only (Ns ignored)."""

    seq = seq.upper()

    valid = [b for b in seq if b in VALID]

    if not valid:

        return float("nan")

    gc = sum(1 for b in valid if b in {"G", "C"})

    return 100.0 * gc / len(valid)



def normalize_species(sp: str) -> str:

    """

    Normalize/merge known aliases:

      - pfoki / Pfokini -> pfoki

      - multimicronucleatum / Pmulti -> multimicronucleatum

      - pquin / Pquin -> pquin

    Also extracts canonical species tokens if present.

    """

    s = sp.strip().lower().replace("-", "_").replace(".", "_")



    # alias merges

    if "pfok" in s:

        return "pfoki"

    if "multimicronucleatum" in s or "pmulti" in s:

        return "multimicronucleatum"

    if "pquin" in s:

        return "pquin"



    # keep canonical token if present anywhere

    for tok in [

        "biaurelia","bursaria","caudatum","tetraurelia","triaurelia","primaurelia",

        "sexaurelia",

        "quadecaurelia","decaurelia","dodecaurelia","octaurelia","pentaurelia",

        "tredecaurelia","novaurelia","sonneborni","jenningsi","gigas","tetrahymena",

        "pcalkinski","psche"

    ]:

        if tok in s:

            return tok



    return s



def detect_species(sample_name: str) -> str:

    """

    Species can appear as:

      - sample__species  (preferred)

      - species_...

      - ..._species

    We search for known tokens anywhere if '__' not present.

    """

    name = sample_name.strip().strip("/")



    # 1) Prefer suffix after "__"

    if "__" in name:

        return normalize_species(name.split("__", 1)[1])



    low = name.lower().replace("-", "_").replace(".", "_")



    # 2) Alias shortcuts

    if "pfok" in low:

        return "pfoki"

    if "multimicronucleatum" in low or "pmulti" in low:

        return "multimicronucleatum"

    if "pquin" in low:

        return "pquin"



    # 3) Search known tokens anywhere in the string

    tokens = [

        "biaurelia","bursaria","caudatum","tetraurelia","triaurelia","primaurelia",

        "sexaurelia",

        "quadecaurelia","decaurelia","dodecaurelia","octaurelia","pentaurelia",

        "tredecaurelia","novaurelia","sonneborni","jenningsi","gigas","tetrahymena",

        "pcalkinski","psche","pquin"

    ]

    for t in tokens:

        if t in low:

            return normalize_species(t)



    # 4) Fallback: first token

    return normalize_species(low.split("_")[0] if "_" in low else low)



# ---------- per-sample GC (mean across genes; each gene equal weight) ----------

sample_mean_gc = {}   # sample -> mean GC% across genes

sample_species = {}   # sample -> species

sample_n_genes = {}   # sample -> number of gene files used



for sample_dir in sorted(p for p in ROOT.iterdir() if p.is_dir()):

    sample = sample_dir.name

    species = detect_species(sample)

    gcs = []



    for fna in sorted(sample_dir.glob("*.fna")):

        seq = ""

        with open(fna) as fh:

            for line in fh:

                if not line.startswith(">"):

                    seq += line.strip()



        gcp = gc_percent(seq)

        if gcp == gcp:  # not NaN

            gcs.append(gcp)



    if gcs:

        sample_mean_gc[sample] = mean(gcs)

        sample_species[sample] = species

        sample_n_genes[sample] = len(gcs)



with open(OUT_SAMPLE, "w", newline="") as out:

    w = csv.writer(out)

    w.writerow(["sample", "species", "n_genes_used", "sample_mean_gc_percent"])

    for sample in sorted(sample_mean_gc.keys()):

        w.writerow([

            sample,

            sample_species[sample],

            sample_n_genes[sample],

            round(sample_mean_gc[sample], 6)

        ])



print(f"✅ wrote {OUT_SAMPLE}")



# ---------- per-species GC (mean of sample means; each sample equal weight) ----------

by_species = defaultdict(list)



for sample, gc_mean in sample_mean_gc.items():

    by_species[sample_species[sample]].append(gc_mean)



with open(OUT_SPECIES, "w", newline="") as out:

    w = csv.writer(out)

    w.writerow(["species", "n_samples", "species_mean_gc_percent", "min_sample_gc", "max_sample_gc"])

    for sp in sorted(by_species.keys()):

        vals = by_species[sp]

        w.writerow([

            sp,

            len(vals),

            round(mean(vals), 6),

            round(min(vals), 6),

            round(max(vals), 6),

        ])



print(f"✅ wrote {OUT_SPECIES}")


