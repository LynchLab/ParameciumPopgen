#!/usr/bin/env python3
import re
from pathlib import Path

# ---- run from: .../burs_blast/rev_blastp ----
WD = Path.cwd()

# inputs in THIS folder
IDS_PATH = WD / "bursaria_ids.csv"

# gene-level FASTAs in ../diploid_out
DIP = WD.parent / "diploid_out"
CDS_FA  = DIP / "pburs_diploid.genelevel_cds.fa"
PROT_FA = DIP / "pburs_diploid.genelevel_proteins_code6.fa"

# outputs in ../diploid_out/{nuc_seq,prot_seq}
OUT_NUC  = DIP / "nuc_seq"
OUT_PROT = DIP / "prot_seq"
OUT_NUC.mkdir(parents=True, exist_ok=True)
OUT_PROT.mkdir(parents=True, exist_ok=True)


def read_ids(path: Path):
    ids = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            # tolerate csv/tsv/space
            s = re.split(r"[,\t\s]+", s)[0]
            ids.append(s)
    return ids


def p_to_t(pid: str) -> str:
    # PBUR... .P00120006 -> .T00120006
    return re.sub(r"\.P(\d+)\b", r".T\1", pid)


def iter_fasta(path: Path):
    header = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                if header is None:
                    continue
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def chosen_tx_from_header(h: str):
    m = re.search(r"\bchosen_tx=([A-Za-z0-9_.-]+)", h)
    return m.group(1) if m else None


def wrap(seq: str, w: int = 60) -> str:
    return "\n".join(seq[i:i+w] for i in range(0, len(seq), w))


def build_map_by_tx(fa: Path):
    d = {}
    for hdr, seq in iter_fasta(fa):
        tx = chosen_tx_from_header(hdr)
        if tx:
            d[tx] = (hdr, seq)
    return d


def main():
    if not IDS_PATH.exists():
        raise SystemExit(f"Missing {IDS_PATH}")

    if not CDS_FA.exists():
        raise SystemExit(f"Missing {CDS_FA}")

    if not PROT_FA.exists():
        raise SystemExit(f"Missing {PROT_FA}")

    p_ids = read_ids(IDS_PATH)
    t_ids = [p_to_t(x) for x in p_ids]

    cds_by_tx = build_map_by_tx(CDS_FA)
    prot_by_tx = build_map_by_tx(PROT_FA)

    missing_nuc = []
    missing_prot = []
    wrote_nuc = 0
    wrote_prot = 0

    for pid, tid in zip(p_ids, t_ids):
        # nucleotide
        if tid in cds_by_tx:
            hdr, seq = cds_by_tx[tid]
            out = OUT_NUC / f"{tid}.fa"
            with out.open("w") as w:
                w.write(f">{tid} source_id={pid} {hdr}\n")
                w.write(wrap(seq) + "\n")
            wrote_nuc += 1
        else:
            missing_nuc.append((pid, tid))

        # protein
        if tid in prot_by_tx:
            hdr, seq = prot_by_tx[tid]
            out = OUT_PROT / f"{tid}.fa"
            with out.open("w") as w:
                w.write(f">{tid} source_id={pid} {hdr}\n")
                w.write(wrap(seq) + "\n")
            wrote_prot += 1
        else:
            missing_prot.append((pid, tid))

    print(f"Run dir: {WD}")
    print(f"IDs read: {len(p_ids)}")
    print(f"Nuc written:  {wrote_nuc} -> {OUT_NUC}")
    print(f"Prot written: {wrote_prot} -> {OUT_PROT}")

    if missing_nuc:
        mp = WD / "missing_nuc.txt"
        with mp.open("w") as w:
            for pid, tid in missing_nuc:
                w.write(f"{pid}\t{tid}\n")
        print(f"Missing nuc:  {len(missing_nuc)} (see {mp})")

    if missing_prot:
        mp = WD / "missing_prot.txt"
        with mp.open("w") as w:
            for pid, tid in missing_prot:
                w.write(f"{pid}\t{tid}\n")
        print(f"Missing prot: {len(missing_prot)} (see {mp})")


if __name__ == "__main__":
    main()

