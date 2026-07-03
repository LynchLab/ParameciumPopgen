#!/usr/bin/env python3

import csv

inp = "caud_grouped_by_ref_tetr.csv"
out = "caud_ids_with_tetr.csv"

seen = set()

with open(inp, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        ref = (row.get("ref_tetr_orthologs") or "").strip()
        cauds = (row.get("caud_proteins") or "").strip()

        # keep only rows where tetra orthologs exist AND caud list exists
        if not ref or not cauds:
            continue

        for c in cauds.split(","):
            c = c.strip()
            if c:
                seen.add(c)

with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["caud_protein"])
    for c in sorted(seen):
        w.writerow([c])

print(f"✅ wrote {out} ({len(seen)} unique caud IDs)")

