#!/usr/bin/env python3

import os

import multiprocessing

from pathlib import Path

from Bio import AlignIO

from Bio.SeqRecord import SeqRecord

from Bio.Seq import Seq

from Bio.Data import CodonTable

from collections import Counter



# === ROOT WITH ALL SPECIES ===

PNPS_ROOT = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pn_ps_new")



# === Genetic Code Setup (Paramecium mitochondria) ===

genetic_code_id = 4  # Paramecium mitochondria

ctab = CodonTable.unambiguous_dna_by_id[genetic_code_id]

start_codons = [c for c in ctab.start_codons if c != 'TTA'] + ['GTA']

stop_codons = ctab.stop_codons



print(f"🧬 Using Genetic Code ID: {genetic_code_id}")

print(f"✅ Start codons: {start_codons}")

print(f"🚫 Stop codons: {stop_codons}")



# === Consensus Helper (same logic as your original) ===

def consensus(aln, minc=0.6):

    ncol = aln.get_alignment_length()

    cons = ''

    for i in range(ncol):

        column = [record.seq[i] for record in aln]

        freq = Counter(column)

        top, count = freq.most_common(1)[0]

        if count / len(column) >= minc:

            cons += top

        else:

            cons += 'N'

    return SeqRecord(Seq(cons), id="Consensus")



# === Codon Splitter (same) ===

def n2c(nuc_seq):

    seq = str(nuc_seq).upper()

    return [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]



# === Codon Degeneracy Table Builder (same logic) ===

def sitedegen(ct=ctab, nucs='ACGT'):

    ixs = range(3)

    table = ct.forward_table

    codons = table.keys()

    out = {k: () for k in codons}

    for codon in codons:

        aa = table[codon]

        f = [0 for _ in ixs]

        for i in ixs:

            altc = [''.join([n if x == i else codon[x] for x in ixs]) for n in nucs]

            alta = [table[c] for c in altc if c not in stop_codons and c in table]

            # same as your code: size of the most-common synonymous AA set

            f[i] = Counter(alta).most_common(1)[0][1] if alta else 0

        out[codon] = f

    return out



degen_tab = sitedegen()

print("🧬 Codon degeneracy table built.")



# === Compute Degeneracy for a Single Alignment (same logic) ===

def get_dseq(fin, fout, dtab=degen_tab):

    if os.path.exists(fout):

        print(f"⏩ Skipping (already exists): {fout}")

        return

    try:

        aln = AlignIO.read(fin, 'fasta')

    except ValueError:

        print(f"❌ Could not parse alignment: {fin}")

        return



    conseq = consensus(aln, minc=0.6)

    codons = n2c(conseq.seq)

    dlist = [dtab.get(i, ['0', '0', '0']) for i in codons]

    dseq = ''.join(str(i) for triplet in dlist for i in triplet)



    with open(fout, 'w') as f:

        f.write(dseq + '\n')

    print(f"✅ Wrote degeneracy: {fout}")



# === Helper to derive a gene name from filename (non-invasive) ===

def gene_from_fname(fname: str) -> str:

    if fname.endswith("_codon_aligned.fna"):

        return fname[:-len("_codon_aligned.fna")]

    if fname.endswith(".fna"):

        return fname[:-4]

    return Path(fname).stem



def process_species(species_dir: Path):

    codon_dir = species_dir / "codon_alignment"

    if not codon_dir.is_dir():

        print(f"[SKIP] {species_dir.name}: no codon_alignment/")

        return



    # Outputs go directly here (no alignment_input subfolder)

    output_dir = species_dir / "degeneracy_result"

    output_dir.mkdir(parents=True, exist_ok=True)



    # Collect inputs from codon_alignment

    alignment_files = [f for f in os.listdir(codon_dir) if f.endswith(".fna")]

    if not alignment_files:

        print(f"[SKIP] {species_dir.name}: no .fna in codon_alignment/")

        return



    fins  = [str(codon_dir / f) for f in alignment_files]

    fouts = [str(output_dir / f"{gene_from_fname(f)}.txt") for f in alignment_files]

    flist = list(zip(fins, fouts))



    print(f"🚀 Starting degeneracy analysis for {species_dir.name}: {len(flist)} files")

    with multiprocessing.Pool(processes=6) as pool:

        pool.starmap(get_dseq, flist)



def main():

    for species_dir in sorted(p for p in PNPS_ROOT.iterdir() if p.is_dir()):

        process_species(species_dir)

    print("🎉 All done.")



if __name__ == "__main__":

    main()


