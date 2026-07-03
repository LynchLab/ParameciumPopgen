import math, csv

from pathlib import Path



# ✅ Use the actual folder you're in

root = Path(".")

tsv = root / "m0_m3_lnL.tsv"

csv_out = root / "m0_m3_lrt_results.csv"



rows = []

with tsv.open() as f:

    header = f.readline()

    for line in f:

        gene, l0, l1, stat = line.strip().split(" ")

        if gene == "gene":  # skip header if present again

            continue

        try:

            x = float(stat)

            # clamp negatives to 0 for the LRT (standard practice)

            if x < 0:

                x = 0.0

            # Mixture ½·χ²(df=1) p-value

            p_mix = 0.5 * math.erfc((x / 2.0) ** 0.5)

            sig5 = "YES" if x >= 2.71 else "NO"

            sig1 = "YES" if x >= 5.41 else "NO"

            rows.append([gene, l0, l1, f"{x:.6f}", f"{p_mix:.6g}", sig5, sig1])

        except ValueError:

            rows.append([gene, l0, l1, "NA", "NA", "NA", "NA"])



with csv_out.open("w", newline="") as f:

    w = csv.writer(f)

    w.writerow([

        "gene",

        "lnL_M8a",

        "lnL_M8",

        "stat_2deltaLL",

        "p_mixture_half_chi2_df1",

        "sig_5pct",

        "sig_1pct"

    ])

    w.writerows(rows)



print(f"✅ LRT results written to: {csv_out.resolve()}")

