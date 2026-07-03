#!/usr/bin/env python3
import os
import re
import csv
import glob
from typing import Optional, Tuple

# Example line in .mlc:
# p0 =   0.99885  p =   0.46696 q =   9.96364
PQ_RE = re.compile(r"\bp\s*=\s*([0-9]*\.?[0-9]+)\s+q\s*=\s*([0-9]*\.?[0-9]+)\b")

def extract_pq_from_mlc(mlc_path: str) -> Optional[Tuple[float, float]]:
    """
    Return (p, q) from an M8 .mlc file. Looks for the first occurrence of 'p = ... q = ...'.
    """
    with open(mlc_path, "r", errors="replace") as f:
        for line in f:
            m = PQ_RE.search(line)
            if m:
                p = float(m.group(1))
                q = float(m.group(2))
                return p, q
    return None

def main():
    base_dir = os.getcwd()  # run from M8_vs_M8a
    out_csv = os.path.join(base_dir, "m8_pq_omega_new.csv")

    rows = []
    skipped = []

    # Looks like: OMA01173/M8/OMA01173_M8.mlc
    mlc_paths = sorted(glob.glob(os.path.join(base_dir, "OMA*/M8/*_M8.mlc")))

    for mlc in mlc_paths:
        # gene folder is OMAxxxxx
        parts = mlc.split(os.sep)
        # find the first component that starts with "OMA"
        gene = next((p for p in parts if p.startswith("OMA")), None)
        if gene is None:
            skipped.append((mlc, "Could not infer gene (OMA folder)"))
            continue

        pq = extract_pq_from_mlc(mlc)
        if pq is None:
            skipped.append((gene, "Could not find p/q in .mlc"))
            continue

        p, q = pq
        omega = p / (p + q) if (p + q) != 0 else float("nan")

        rows.append({
            "gene": gene,
            "p": p,
            "q": q,
            "omega": omega
        })

    # write output
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "p", "q", "omega"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # also write a small skipped log for debugging
    if skipped:
        skipped_path = os.path.join(base_dir, "m8_pq_omega.skipped.tsv")
        with open(skipped_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["item", "reason"])
            for item, reason in skipped:
                w.writerow([item, reason])

    print(f"Wrote: {out_csv}  (rows={len(rows)})")
    if skipped:
        print(f"Skipped: {len(skipped)}  -> m8_pq_omega.skipped.tsv")

if __name__ == "__main__":
    main()

