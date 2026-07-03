#!/usr/bin/env python3

import sys



pcaud_file = "pcaud_prots_w_tetrahymena_etc_homologs.txt"

hog_file   = "root_hogs_protein_ids.tsv"

out_file   = "pcaud_paralogs_orthologs.tsv"



# read PCAU protein IDs

with open(pcaud_file) as f:

    pcaud_ids = {line.strip() for line in f if line.strip()}



results = {pid: [] for pid in pcaud_ids}



with open(hog_file) as f:

    for line in f:

        cols = line.rstrip().split("\t")

        hog  = cols[0]

        prots = cols[1:]



        for p in prots:

            # extract raw protein ID before "__[species]"

            pid = p.split("__")[0]

            if pid in pcaud_ids:

                # collect all other proteins in same HOG

                others = [x.split("__")[0] for x in prots if x.split("__")[0] != pid]

                results[pid].extend(others)



# write output

with open(out_file, "w") as out:

    out.write("PCAU_protein\torthologs_paralogs\n")

    for pid, hits in results.items():

        uniq_hits = sorted(set(hits))

        out.write(f"{pid}\t{','.join(uniq_hits)}\n")



print(f"✅ Written: {out_file}")


