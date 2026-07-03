import pandas as pd

import numpy as np

from pathlib import Path



# folder containing your CSVs

folder = Path(".")  # change if needed



rows = []

for fp in sorted(folder.glob("average_genetic_distance_4fold_*.csv")):

    species = fp.stem.replace("average_genetic_distance_4fold_", "")

    df = pd.read_csv(fp, index_col=0)



    # Ensure numeric

    df = df.apply(pd.to_numeric, errors="coerce")



    # Use only the upper triangle (i<j), exclude self and zeros/NaN

    iu = np.triu_indices_from(df, k=1)

    vals = df.values[iu]

    vals = vals[~np.isnan(vals)]

    vals = vals[vals > 0]



    mean_ds = float(np.mean(vals)) if vals.size else np.nan

    n_pairs = int(vals.size)



    rows.append({"species": species, "mean_dS": mean_ds, "n_pairs": n_pairs})



out = pd.DataFrame(rows).sort_values("species")

print(out)



# save summary

out.to_csv("average_dS_per_species.csv", index=False)

print("\nSaved: average_dS_per_species.csv")


