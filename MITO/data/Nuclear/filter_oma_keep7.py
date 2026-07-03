#!/usr/bin/env python3

KEEP_SPECIES = {
    "biau",
    "burs_hap",
    "caud",
    "mult",
    "prim_Ir4",
    "tetr",
    "sexa"
}

input_file = "OMA_rows_for_caud_with_tetr_onlyPBURShap.txt"
output_file = "OMA_rows_keep7.txt"

with open(input_file) as fin, open(output_file, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue

        fields = line.split()
        oma_id = fields[0]

        kept = [oma_id]

        for field in fields[1:]:
            # species is everything before the first colon
            species = field.split(":", 1)[0]
            if species in KEEP_SPECIES:
                kept.append(field)

        fout.write("\t".join(kept) + "\n")

print(f"✅ Wrote filtered file: {output_file}")

