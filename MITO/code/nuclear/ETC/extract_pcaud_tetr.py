#!/usr/bin/env python3

import csv

input_tsv = "pcaud_paralogs_orthologs.tsv"
output_csv = "caud_ref_tetr_only.csv"

with open(input_tsv, newline="") as fin, open(output_csv, "w", newline="") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    writer = csv.writer(fout)
    
    # Write header
    writer.writerow(["caud_protein", "ref_tetr_orthologs"])
    
    for row in reader:
        caud = row["PCAU_protein"]
        
        orthos = row["orthologs_paralogs"].split(",")
        ref_tetr = [x for x in orthos if x.startswith("PTET.")]
        
        if ref_tetr:
            writer.writerow([caud, ",".join(ref_tetr)])

print("✅ Wrote caud_ref_tetr_only.csv")

