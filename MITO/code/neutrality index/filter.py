#!/usr/bin/env python3

import pandas as pd



IN  = "combined_normalized_internal_nodes_summary.csv"

OUT = "ancestral_dnds_highprob_filterallsites_0.60.csv"



keep = [

    "gene","species","node","dn","ds",

    "synonymous_sites","non_synonymous_sites",

    "Normalized_dN","Normalized_dS","Normalized_dN/dS",

]



df = pd.read_csv(IN)

# clean up any stray spaces in headers from wrapping

df.columns = [c.strip() for c in df.columns]



missing = [c for c in keep if c not in df.columns]

if missing:

    raise SystemExit(f"Missing expected columns: {missing}")



df[keep].to_csv(OUT, index=False)

print(f"✅ wrote {OUT} with {len(df)} rows and {len(keep)} columns")


