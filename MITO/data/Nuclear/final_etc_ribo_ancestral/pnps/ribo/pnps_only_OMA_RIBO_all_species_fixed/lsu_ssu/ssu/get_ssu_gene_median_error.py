#!/usr/bin/env python3

import pandas as pd
import numpy as np

INFILE = "SSU_nuclear_gene_median_pnps.csv"
OUTFILE = "SSU_nuclear_gene_median_pnps_errorbars.csv"

df = pd.read_csv(INFILE)
df.columns = [c.strip() for c in df.columns]

df["median_pN_pS"] = pd.to_numeric(df["median_pN_pS"], errors="coerce")
df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["Gene", "median_pN_pS"])

vals = df["median_pN_pS"].dropna()

median = vals.median()
q1 = vals.quantile(0.25)
q3 = vals.quantile(0.75)

error_low = median - q1
error_high = q3 - median

out = pd.DataFrame([{
    "group": "Nuclear_SSU_rpS",
    "n_genes": df["Gene"].nunique(),
    "n_values": len(vals),
    "median_pN_pS": median,
    "q1": q1,
    "q3": q3,
    "error_low": error_low,
    "error_high": error_high,
    "min": vals.min(),
    "max": vals.max()
}])

out.to_csv(OUTFILE, index=False)

print(f"Wrote: {OUTFILE}")
print(out.to_string(index=False))
print()
print(f"Formatted: {median:.4f} -{error_low:.4f} / +{error_high:.4f}")
