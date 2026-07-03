#!/usr/bin/env python3

from pathlib import Path

from collections import Counter

from Bio import AlignIO



VALID = set("ATGC")



EXCLUDE = {

    "gc_per_species_all.csv",

    "gc_per_species.csv",

    "gc_per_species_intergenic.py",

    "nuc_div.py",

    "pi_results",

    "mafft.sh",

    "scripts",

}



def pi_for_column(bases):

    m = len(bases)

    if m < 2:

        return 0.0

    c = Counter(bases)

    pairs = m * (m - 1) / 2

    same = sum(v * (v - 1) / 2 for v in c.values())

    return (pairs - same) / pairs



def compute_alignment_pi(aln_file, out_tsv):

    try:

        aln = AlignIO.read(aln_file, "fasta")

    except Exception:

        return False, 0, 0.0



    # 🔴 Ignore regions with <2 samples

    if len(aln) < 2:

        return False, 0, 0.0



    L = aln.get_alignment_length()

    total_pi = 0.0

    n_sites = 0



    with open(out_tsv, "w") as out:

        out.write("Position\tPi\n")

        for i in range(L):

            col = [rec.seq[i].upper() for rec in aln]

            bases = [b for b in col if b in VALID]

            pi = pi_for_column(bases)

            out.write(f"{i+1}\t{pi}\n")

            total_pi += pi

            n_sites += 1



    return True, n_sites, total_pi



def main():

    root = Path(".")



    print("species\tregions_used\ttotal_sites\tsite_weighted_mean_pi")



    for species_dir in sorted(p for p in root.iterdir() if p.is_dir() and p.name not in EXCLUDE):

        mafft_dir = species_dir / "mafft"

        if not mafft_dir.exists():

            continue



        pi_dir = species_dir / "pi"

        pi_dir.mkdir(exist_ok=True)



        total_sites = 0

        total_pi = 0.0

        regions_used = 0



        for aln in sorted(mafft_dir.glob("*.mafft.fa")):

            region = aln.stem.replace(".mafft", "")

            out_tsv = pi_dir / f"{region}_pi_per_site.tsv"



            ok, n_sites, sum_pi = compute_alignment_pi(aln, out_tsv)

            if ok:

                regions_used += 1

                total_sites += n_sites

                total_pi += sum_pi

            else:

                if out_tsv.exists():

                    out_tsv.unlink()



        mean_pi = total_pi / total_sites if total_sites else float("nan")



        # per-species summary

        with open(pi_dir / "species_pi_summary.txt", "w") as out:

            out.write(f"species\t{species_dir.name}\n")

            out.write(f"regions_used\t{regions_used}\n")

            out.write(f"total_sites\t{total_sites}\n")

            out.write(f"site_weighted_mean_pi\t{mean_pi}\n")



        print(f"{species_dir.name}\t{regions_used}\t{total_sites}\t{mean_pi}")



if __name__ == "__main__":

    main()


