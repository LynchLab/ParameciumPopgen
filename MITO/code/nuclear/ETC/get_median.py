import pandas as pd
import numpy as np
from pathlib import Path

INFILE = Path("M8_beta_pq_omega.csv")
OUTFILE = Path("M8_beta_pq_omega_median_error.csv")

df = pd.read_csv(INFILE)

# Convert numeric columns
numeric_cols = ["p", "q", "omega_w_gt_1", "omega"]

for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.replace([np.inf, -np.inf], np.nan)

rows = []

for col in numeric_cols:
    vals = df[col].dropna()

    median = vals.median()
    q1 = vals.quantile(0.25)
    q3 = vals.quantile(0.75)

    rows.append({
        "metric": col,
        "n_values": len(vals),
        "median": median,
        "q1": q1,
        "q3": q3,
        "error_low": median - q1,
        "error_high": q3 - median,
        "min": vals.min(),
        "max": vals.max()
    })

out = pd.DataFrame(rows)
out.to_csv(OUTFILE, index=False)

print(f"wrote {OUTFILE}")
print(out.to_string(index=False))
