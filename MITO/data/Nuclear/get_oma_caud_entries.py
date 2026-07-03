#!/usr/bin/env python3

import csv

caud_list_file = "caud_ids_with_tetr.csv"
oma_file = "OrthologousGroups.txt"
out_file = "OMA_rows_for_caud_with_tetr.txt"

# read caud IDs
caud_ids = set()
with open(caud_list_file, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        caud_ids.add(row["caud_protein"].strip())

print(f"Loaded {len(caud_ids)} caudatum IDs")

with open(oma_file) as fin, open(out_file, "w") as fout:
    for line in fin:
        # skip comments
        if line.startswith("#"):
            continue

        # check if any caud ID is present in this OMA row
        if any(f"caud:{cid}" in line for cid in caud_ids):
            fout.write(line)

print(f"✅ Wrote {out_file}")

