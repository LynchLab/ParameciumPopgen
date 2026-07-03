#!/usr/bin/env python3
# filter_oma_by_bursaria_7species_with_missing.py

bursaria_file = "bursaria_list.csv"
oma_file = "OrthologousGroups.txt"
out_file = "OrthologousGroups_with_bursaria_7species.txt"
missing_file = "bursaria_ids_not_found.txt"

keep_species = ["burs", "prim_Ir4", "tetr", "sexa", "caud", "biau", "mult"]

# read PBUR IDs
with open(bursaria_file) as f:
    burs_ids = {line.strip() for line in f if line.strip()}

found_burs_ids = set()

matched = 0
written = 0

with open(oma_file) as fin, open(out_file, "w") as fout:
    for line in fin:
        if line.startswith("#") or not line.strip():
            continue

        # only keep rows containing any PBUR.Dd1.1.Pxxxxx from bursaria_list.csv
        row_burs_hits = [bid for bid in burs_ids if bid in line]
        if not row_burs_hits:
            continue

        matched += 1
        found_burs_ids.update(row_burs_hits)

        parts = line.rstrip("\n").split("\t")
        oma_id = parts[0]
        members = parts[1:]

        # collect only desired species entries
        keep_map = {}
        for m in members:
            if ":" not in m:
                continue
            sp, rest = m.split(":", 1)
            if sp in keep_species:
                keep_map[sp] = m

        # enforce burs presence
        if "burs" not in keep_map:
            continue

        out_parts = [oma_id] + [keep_map[sp] for sp in keep_species if sp in keep_map]
        fout.write("\t".join(out_parts) + "\n")
        written += 1

# write missing PBUR IDs
missing_ids = sorted(burs_ids - found_burs_ids)
with open(missing_file, "w") as f:
    for bid in missing_ids:
        f.write(bid + "\n")

print(f"Matched rows (by PBUR list): {matched}")
print(f"Wrote rows (7-species subset): {written} -> {out_file}")
print(f"PBUR IDs not found in any OMA group: {len(missing_ids)} -> {missing_file}")

