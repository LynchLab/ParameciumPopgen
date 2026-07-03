#!/usr/bin/env python3
import os
from collections import defaultdict

# -------------------------------------------------
# RUN FROM:
# /data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/rev_blastp
# -------------------------------------------------

BASE_DIR = os.getcwd()

BURS_IDS = os.path.join(BASE_DIR, "bursaria_ids.csv")
HOGS_TSV = os.path.join(BASE_DIR, "..", "..", "etc", "root_hogs_protein_ids.tsv")
OUT_TSV  = os.path.join(BASE_DIR, "burs_to_caud.tsv")


def read_burs_ids(path):
    ids = []
    seen = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.lower() in {"pburs_protein_id", "bursaria_id", "id", "qseqid"}:
                continue
            tok = s.split()[0]
            if tok not in seen:
                seen.add(tok)
                ids.append(tok)
    return ids


def main():
    if not os.path.exists(BURS_IDS):
        raise FileNotFoundError(f"Missing: {BURS_IDS}")
    if not os.path.exists(HOGS_TSV):
        raise FileNotFoundError(f"Missing: {HOGS_TSV}")

    burs_ids = read_burs_ids(BURS_IDS)

    # PBUR -> set(PCAU)
    burs_to_caud = defaultdict(set)

    with open(HOGS_TSV) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                continue

            pburs = []
            pcauds = []

            for tok in fields[1:]:
                if not tok:
                    continue

                core = tok.split("__", 1)[0]

                # bursaria (ignore burs_hap)
                if core.startswith("PBUR.") and "__[burs]" in tok:
                    pburs.append(core)

                # caudatum
                if core.startswith("PCAU.") and "__[caud]" in tok:
                    pcauds.append(core)

            if not pburs:
                continue

            for pb in pburs:
                for pc in pcauds:
                    burs_to_caud[pb].add(pc)

    # -------------------------
    # WRITE OUTPUT
    # -------------------------
    with open(OUT_TSV, "w") as out:
        out.write("PBUR_ID\tPCAU_IDs\n")

        for pb in burs_ids:
            cauds = sorted(burs_to_caud.get(pb, []))
            if cauds:
                out.write(f"{pb}\t{','.join(cauds)}\n")
            else:
                out.write(f"{pb}\tNA\n")

    mapped = sum(1 for pb in burs_ids if pb in burs_to_caud and burs_to_caud[pb])
    print(f"Mapped {mapped}/{len(burs_ids)} PBUR IDs")
    print(f"Wrote: {OUT_TSV}")


if __name__ == "__main__":
    main()

