#!/usr/bin/env python3

import os

import pandas as pd

from pathlib import Path

from Bio import AlignIO

from Bio.Align import MultipleSeqAlignment

from Bio.SeqRecord import SeqRecord

from Bio.Seq import Seq

from utils.custom_classes import MSA



# Root with all species

PNPS_ROOT = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pn_ps_new")



# Parameters (match your MSA usage)

alignment_ext = "_codon_aligned.fna"

cex = "N-"   # exclude chars

t = 0.1      # threshold for MSA.pdistmat



def get_gene_name(fname: str) -> str:

    return fname[:-len(alignment_ext)] if fname.endswith(alignment_ext) else Path(fname).stem



def extract_4fold_sites(aln_path: Path, degen_path: Path, out_dir: Path, gene_name: str):

    try:

        aln = AlignIO.read(str(aln_path), "fasta")

    except Exception as e:

        print(f"❌ read alignment failed: {aln_path} — {e}")

        return None, None

    try:

        degen = Path(degen_path).read_text().strip()

    except Exception as e:

        print(f"❌ read degeneracy failed: {degen_path} — {e}")

        return None, None



    target_char = '4' if '4' in degen else '3'

    if target_char == '3' and '4' not in degen:

        print(f"[INFO] {gene_name}: degeneracy appears 0–3 (synonymous-count). Using '3' as 4-fold sites.")

    keep_indices = [i for i, d in enumerate(degen) if d == target_char]

    if not keep_indices:

        print(f"[WARN] {gene_name}: no 4-fold indices found in {degen_path}")

        return None, None



    L = min(aln.get_alignment_length(), len(degen))

    keep_indices = [i for i in keep_indices if i < L]



    filtered = []

    for rec in aln:

        s = str(rec.seq)[:L]

        filtered_seq = ''.join(s[i] for i in keep_indices)

        filtered.append(SeqRecord(Seq(filtered_seq), id=rec.id, description=""))



    out_dir.mkdir(parents=True, exist_ok=True)

    filt_fa = out_dir / f"{gene_name}_4fold_filtered.fna"

    AlignIO.write(MultipleSeqAlignment(filtered), str(filt_fa), "fasta")



    idx_path = out_dir / f"{gene_name}_4fold_indices.txt"

    idx_path.write_text(",".join(map(str, keep_indices)))

    print(f"✅ {gene_name}: wrote 4-fold filtered alignment → {filt_fa}")

    return filt_fa, idx_path



def write_distance_matrices(filtered_fa: Path, out_dir: Path, gene_name: str):

    msa = MSA(str(filtered_fa), sfx=".fna")

    dmat = msa.pdistmat(cex, t)



    full_csv = out_dir / f"{gene_name}_distance_matrix.csv"

    dmat.to_csv(full_csv)

    print(f"📏 {gene_name}: full distance matrix → {full_csv}")



    labels = msa.seqids

    n = len(labels)

    lower = pd.DataFrame(index=labels, columns=labels, dtype=object)

    for i in range(n):

        for j in range(i + 1):

            val = dmat.iloc[i, j]

            if pd.isna(val):

                val = dmat.iloc[j, i]

            lower.iloc[i, j] = round(float(val), 5)

    lower_csv = out_dir / f"{gene_name}_lower_triangle_labeled.csv"

    lower.to_csv(lower_csv)

    print(f"📁 {gene_name}: lower triangle → {lower_csv}")



def write_species_average(out_root: Path):

    """Aggregate all *_lower_triangle_labeled.csv under out_root/<gene>/ into a single average matrix."""

    pairwise = {}  # (s1,s2) -> [vals]

    genes = [p for p in out_root.iterdir() if p.is_dir()]

    for gdir in genes:

        gene = gdir.name

        f = gdir / f"{gene}_lower_triangle_labeled.csv"

        if not f.exists():

            continue

        try:

            df = pd.read_csv(f, index_col=0)

        except Exception as e:

            print(f"[WARN] could not read {f}: {e}")

            continue

        samples = list(df.index)

        for i, s1 in enumerate(samples):

            for j in range(i + 1):

                s2 = samples[j]

                val = df.iloc[i, j]

                if pd.isna(val):

                    continue

                key = tuple(sorted((s1, s2)))

                pairwise.setdefault(key, []).append(float(val))



    if not pairwise:

        print(f"[SKIP] no pairwise values found in {out_root}")

        return



    all_samples = sorted({s for pair in pairwise for s in pair})

    avg_mat = pd.DataFrame(index=all_samples, columns=all_samples, dtype=float)

    for s in all_samples:

        avg_mat.at[s, s] = 0.0

    for (s1, s2), vals in pairwise.items():

        avg = sum(vals) / len(vals)

        avg_mat.at[s1, s2] = avg

        avg_mat.at[s2, s1] = avg



    out_csv = out_root / "average_genetic_distance_4fold_pergene.csv"

    avg_mat.to_csv(out_csv)

    print(f"🧮 wrote species average matrix → {out_csv}")



def main():

    for species_dir in sorted(p for p in PNPS_ROOT.iterdir() if p.is_dir()):

        codon_dir = species_dir / "codon_alignment"

        degen_dir = species_dir / "degeneracy_result"

        out_root  = species_dir / "genetic_distance_4fold"



        if not codon_dir.is_dir():

            print(f"[SKIP] {species_dir.name}: missing codon_alignment/")

            continue

        if not degen_dir.is_dir():

            print(f"[SKIP] {species_dir.name}: missing degeneracy_result/")

            continue



        print(f"\n🔬 Processing species: {species_dir.name}")

        fa_files = sorted(f for f in os.listdir(codon_dir) if f.endswith(alignment_ext))

        if not fa_files:

            print(f"[SKIP] {species_dir.name}: no {alignment_ext} files")

            continue



        for fa in fa_files:

            gene = get_gene_name(fa)

            aln_path   = codon_dir / fa

            degen_path = degen_dir / f"{gene}.txt"

            gene_out   = out_root / gene



            if not degen_path.exists():

                print(f"[WARN] {species_dir.name}:{gene}: missing degeneracy file {degen_path}")

                continue



            filt_fa, _ = extract_4fold_sites(aln_path, degen_path, gene_out, gene)

            if filt_fa is None:

                continue

            try:

                write_distance_matrices(filt_fa, gene_out, gene)

            except Exception as e:

                print(f"⚠️  {species_dir.name}:{gene}: distance calc failed — {e}")



        # After all genes for this species, compute the average matrix

        write_species_average(out_root)



    print("\n🎉 Done: genetic_distance_4fold created per species with per-gene outputs and species-level average.")



if __name__ == "__main__":

    main()


