#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# === Roots (top-level dirs) ===
ROOT = Path(".")
PROT_DIR = ROOT / "mafft_prot"        # input: OMA*.aln.fa  (AA alignments from MAFFT)
NUC_DIR  = ROOT / "nucleotide_seq"    # input: OMA*.fa      (nucleotide CDS, same headers)
OUT_DIR  = ROOT / "codon_alignment"   # output: OMA*.codon.aln.fna

# Include TGA as acceptable terminal codon (you asked: "if TGA put that")
STOP_CODONS = {"TAA", "TAG", "TGA"}

def clean_nuc(seq: str, rec_id: str, gene: str) -> str:
    """Uppercase, T for U, strip gaps/spaces, trim to multiple of 3, ensure terminal stop codon."""
    s = str(seq).upper().replace("U", "T").replace("-", "").replace(" ", "")
    if len(s) < 3:
        return s

    rem = len(s) % 3
    if rem != 0:
        print(f"[WARN] {gene}:{rec_id} length {len(s)} not %3==0 -> trimming {rem} base(s)", file=sys.stderr)
        s = s[: len(s) - rem]

    last = s[-3:] if len(s) >= 3 else ""
    if last not in STOP_CODONS:
        print(f"[INFO] {gene}:{rec_id} missing terminal stop (got '{last or 'NA'}') -> appending TAA", file=sys.stderr)
        s += "TAA"
    return s

def protein_to_codon_aln(prot_aln: str, nuc_seq: str) -> str:
    """
    Project AA alignment onto nucleotide codons:
      - AA gap '-' -> '---'
      - AA letter consumes next codon (excluding final terminal codon)
      - append terminal codon at end (TAA/TAG/TGA)
    """
    codons = [nuc_seq[i:i+3] for i in range(0, len(nuc_seq), 3)]
    if not codons:
        return "---" * (len(prot_aln) + 1)

    out = []
    ci = 0
    for aa in prot_aln:
        if aa == "-":
            out.append("---")
        else:
            if ci < len(codons) - 1:
                out.append(codons[ci])
                ci += 1
            else:
                out.append("---")  # ran out before the terminal codon
    out.append(codons[-1])  # terminal codon
    return "".join(out)

def validate_codon_alignment(path_fna: Path, gene: str):
    bad_len = []
    bad_stop = []
    for rec in SeqIO.parse(str(path_fna), "fasta"):
        seq = str(rec.seq)
        if set(seq) == {"-"}:
            continue
        if (len(seq) % 3) != 0:
            bad_len.append(rec.id)
        last = seq[-3:] if len(seq) >= 3 else ""
        if last not in STOP_CODONS and last != "---":
            bad_stop.append((rec.id, last))

    if bad_len:
        print(f"[CHECK] {gene}: {len(bad_len)} seqs not multiple of 3 (showing up to 10): {', '.join(bad_len[:10])}", file=sys.stderr)
    else:
        print(f"[CHECK] {gene}: all sequences length %3==0", file=sys.stderr)

    if bad_stop:
        preview = ", ".join(f"{sid}:{last}" for sid, last in bad_stop[:10])
        print(f"[CHECK] {gene}: {len(bad_stop)} seqs not ending TAA/TAG/TGA (showing up to 10): {preview}", file=sys.stderr)
    else:
        print(f"[CHECK] {gene}: all sequences end with TAA/TAG/TGA (or all-gap rows)", file=sys.stderr)

def main():
    if not PROT_DIR.is_dir():
        print(f"[ERROR] Missing protein alignment dir: {PROT_DIR}", file=sys.stderr)
        sys.exit(1)
    if not NUC_DIR.is_dir():
        print(f"[ERROR] Missing nucleotide dir: {NUC_DIR}", file=sys.stderr)
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    prot_files = sorted(PROT_DIR.glob("OMA*.aln.fa"))
    if not prot_files:
        print(f"[ERROR] No OMA*.aln.fa found under {PROT_DIR}", file=sys.stderr)
        sys.exit(1)

    for prot_path in prot_files:
        gene = prot_path.stem.replace(".aln", "")  # OMA41918.aln -> OMA41918
        nuc_path = NUC_DIR / f"{gene}.fa"
        out_path = OUT_DIR / f"{gene}.codon.aln.fna"

        if not nuc_path.exists():
            print(f"[WARN] Missing nucleotide file for {gene}: {nuc_path} -> skipping", file=sys.stderr)
            continue

        # Load nuc FASTA into dict by record.id (strip after whitespace)
        nuc_dict = {}
        for rec in SeqIO.parse(str(nuc_path), "fasta"):
            rid = rec.id.split()[0]  # e.g. biau|PBIA... or caud|PCAU...
            nuc_dict[rid] = clean_nuc(rec.seq, rid, gene)

        codon_records = []
        missing = 0

        for prec in SeqIO.parse(str(prot_path), "fasta"):
            pid = prec.id.split()[0]
            prot_aln = str(prec.seq)

            if pid in nuc_dict:
                codon_seq = protein_to_codon_aln(prot_aln, nuc_dict[pid])
            else:
                codon_seq = ("---" * len(prot_aln)) + "---"
                missing += 1

            prec.seq = Seq(codon_seq)
            codon_records.append(prec)

        SeqIO.write(codon_records, str(out_path), "fasta")
        print(f"✅ {gene}: wrote {out_path} (missing nuc for {missing} seqs)")
        validate_codon_alignment(out_path, gene)

    print("All codon alignments completed.")

if __name__ == "__main__":
    main()

