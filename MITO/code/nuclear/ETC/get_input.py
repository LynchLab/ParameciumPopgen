#!/usr/bin/env python3

import re
from pathlib import Path
from typing import Dict, List, Tuple


# ============================================================
# CONFIGURATION
# ============================================================

RAW_BASE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/raw_files"
)

OMA_FILE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/OMA_rows_keep7.txt"
)

OUT_BASE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/dnds/etc/input_dnds"
)

OUT_PROT = OUT_BASE / "protein_seq"
OUT_NUCL = OUT_BASE / "nucleotide_seq"

# If True: only write an OMA fasta if ALL requested species have BOTH protein+nt
# If False: write whatever species have BOTH protein+nt (never protein-only or nt-only)
REQUIRE_ALL = False


# ============================================================
# Species selection rules
# ============================================================

# output_label : OMA tag in the OMA_rows_keep7.txt file
# NOTE: burs_hap IDs are PBUR.* and are stored in raw burs.fa
WANTED = {
    "biau": "biau",
    "burs": "burs_hap",     # fetch PBUR IDs from burs.fa but label as burs
    "caud": "caud",
    "mult": "mult",
    "prim": "prim_Ir4",     # fetch PPRIMP IDs from prim.fa but label as prim
    "sexa": "sexa",
    "tetr": "tetr",
}

# output_label : raw fasta prefix (raw_files/*_seq/<prefix>.fa)
RAW_PREFIX = {
    "biau": "biau",
    "burs": "burs",   # burs_hap PBUR IDs live in burs.fa (as you showed)
    "caud": "caud",
    "mult": "mult",
    "prim": "prim",   # prim_Ir4 IDs live in prim.fa
    "sexa": "sexa",
    "tetr": "tetr",
}


# ============================================================
# FASTA utilities
# ============================================================

def load_fasta_as_dict(fa_path: Path) -> Dict[str, str]:
    """
    Read FASTA into dict: {seq_id: sequence}
    seq_id = first token after '>'
    """
    seqs: Dict[str, str] = {}
    cur_id = None
    cur_seq: List[str] = []

    with fa_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)

        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq)

    return seqs


def write_fasta(out_path: Path, records: List[Tuple[str, str]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


# ============================================================
# OMA parsing
# ============================================================

def parse_oma_row(line: str):
    """
    Returns (oma_id, tag2id) or (None, None) if line is empty/comment.
    """
    line = line.strip()
    if not line or line.startswith("#"):
        return None, None

    parts = re.split(r"\s+", line)
    oma_id = parts[0]

    tag2id: Dict[str, str] = {}
    for tok in parts[1:]:
        if ":" in tok:
            tag, gid = tok.split(":", 1)
            tag2id[tag.strip()] = gid.strip()

    return oma_id, tag2id


# ============================================================
# MAIN
# ============================================================

def main():
    # Load all raw FASTAs once
    prot_db: Dict[str, Dict[str, str]] = {}
    nucl_db: Dict[str, Dict[str, str]] = {}

    for label, prefix in RAW_PREFIX.items():
        prot_fa = RAW_BASE / "protein_seq" / f"{prefix}.fa"
        nucl_fa = RAW_BASE / "nucleotide_seq" / f"{prefix}.fa"

        if not prot_fa.exists():
            raise RuntimeError(f"Missing protein FASTA: {prot_fa}")
        if not nucl_fa.exists():
            raise RuntimeError(f"Missing nucleotide FASTA: {nucl_fa}")

        prot_db[label] = load_fasta_as_dict(prot_fa)
        nucl_db[label] = load_fasta_as_dict(nucl_fa)

    OUT_PROT.mkdir(parents=True, exist_ok=True)
    OUT_NUCL.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    wrote = 0
    skipped_missing_tags = 0
    skipped_require_all = 0

    with OMA_FILE.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            oma_id, tag2id = parse_oma_row(line)
            if oma_id is None:
                continue

            total_rows += 1

            # build records, but only include a species if BOTH protein and nucleotide exist
            prot_records: List[Tuple[str, str]] = []
            nucl_records: List[Tuple[str, str]] = []

            missing_tags = []
            missing_pairs = []

            for out_label, oma_tag in WANTED.items():
                if oma_tag not in tag2id:
                    missing_tags.append(oma_tag)
                    continue

                gene_id = tag2id[oma_tag]
                pseq = prot_db[out_label].get(gene_id)
                nseq = nucl_db[out_label].get(gene_id)

                # Only keep if both exist (this fixes your mismatches)
                if pseq is None or nseq is None:
                    missing_pairs.append(f"{out_label}|{gene_id}")
                    continue

                prot_records.append((f"{out_label}|{gene_id}", pseq))
                nucl_records.append((f"{out_label}|{gene_id}", nseq))

            if REQUIRE_ALL:
                # must have no missing tags and must have all requested species present as pairs
                if missing_tags:
                    skipped_missing_tags += 1
                    continue
                if len(prot_records) != len(WANTED):
                    skipped_require_all += 1
                    continue

            # write even if partial, but guaranteed matched pairs only
            write_fasta(OUT_PROT / f"{oma_id}.fa", prot_records)
            write_fasta(OUT_NUCL / f"{oma_id}.fa", nucl_records)
            wrote += 1

    print("Total OMA rows scanned:", total_rows)
    print("FASTA written:", wrote)
    print("Skipped (missing required tags):", skipped_missing_tags)
    print("Skipped (REQUIRE_ALL failed due to missing protein/nt pair):", skipped_require_all)
    print("Outputs:")
    print("  Protein   ->", OUT_PROT)
    print("  Nucleotide->", OUT_NUCL)
    print("NOTE: Partial outputs only include species with BOTH protein+nucleotide (no mismatches).")


if __name__ == "__main__":
    main()

