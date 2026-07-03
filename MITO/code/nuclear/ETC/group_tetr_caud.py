#!/usr/bin/env python3

import csv
from collections import defaultdict

input_csv = "caud_ref_tetr_only.csv"
output_csv = "caud_grouped_by_ref_tetr.csv"

groups = defaultdict(list)

with open(input_csv, newline="") as fin:
    reader = csv.DictReader(fin)
    for row in reader:
        caud = row["caud_protein"]
        ref_tetr = row["ref_tetr_orthologs"].strip()
        groups[ref_tetr].append(caud)

with open(output_csv, "w", newline="") as fout:
    writer = csv.writer(fout)
    writer.writerow(["ref_tetr_orthologs", "caud_proteins"])
    
    for ref_tetr, cauds in groups.items():
        writer.writerow([ref_tetr, ",".join(sorted(cauds))])

print("✅ Wrote caud_grouped_by_ref_tetr.csv")

