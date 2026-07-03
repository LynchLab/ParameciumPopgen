#!/usr/bin/env python3

from pathlib import Path



BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/recombination"

IN_DIR = Path(BASE) / "aligned_clean"

OUT_DIR = Path(BASE) / "ldhat_seq"

OUT_DIR.mkdir(parents=True, exist_ok=True)



def read_fasta(fp: Path):

    name, buf = None, []

    with fp.open() as f:

        for line in f:

            line = line.strip()

            if not line:

                continue

            if line.startswith(">"):

                if name is not None:

                    yield name, "".join(buf)

                # LDhat recommends <=30 chars

                name = line[1:].split()[0][:30]

                buf = []

            else:

                buf.append(line.replace(" ", ""))

        if name is not None:

            yield name, "".join(buf)



def wrap(s: str, width: int = 2000):

    return "\n".join(s[i:i+width] for i in range(0, len(s), width))



def main():

    files = sorted(IN_DIR.glob("*_aligned.clean.fa"))

    if not files:

        raise SystemExit(f"No *_aligned.clean.fa found in {IN_DIR}")



    for fa in files:

        sp = fa.name.replace("_aligned.clean.fa", "")

        recs = list(read_fasta(fa))

        if not recs:

            print(f"[SKIP] {fa.name} empty")

            continue



        L = len(recs[0][1])

        if any(len(seq) != L for _, seq in recs):

            raise SystemExit(f"[ERROR] Not aligned (unequal lengths): {fa}")



        out_seq = OUT_DIR / f"{sp}.seq"



        with out_seq.open("w") as out:

            out.write(f"{len(recs)} {L} 1\n")  # 1 = haploid/phased

            for name, seq in recs:

                out.write(f">{name}\n{wrap(seq)}\n")



        print(f"[OK] {sp}: wrote {out_seq} (n={len(recs)}, L={L})")



if __name__ == "__main__":

    main()


