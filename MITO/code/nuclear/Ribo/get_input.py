#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ============================================================
# CONFIGURATION
# ============================================================

# Raw multi-FASTA base (USED for all species EXCEPT burs)
RAW_BASE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/raw_files"
)

# OMA rows (7-species subset you created)
OMA_FILE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/OrthologousGroups_with_bursaria_7species.txt"
)

# Bursaria per-gene FASTA folders (USE THESE instead of raw_files for burs)
BURS_GENE_NUCL_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/nuc_seq"
)
BURS_GENE_PROT_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/prot_seq"
)

# Output base
OUT_BASE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/input"
)

OUT_PROT = OUT_BASE / "prot_seq"
OUT_NUCL = OUT_BASE / "nuc_seq"
REPORT_TSV = OUT_BASE / "missing_report.tsv"

# If True: only write an OMA fasta if ALL 7 species are present AND both prot+nucl exist
REQUIRE_ALL_7 = False


# ============================================================
# Species selection rules
# ============================================================

# Final species labels you want in outputs (and in the report)
KEEP_SPECIES = ["burs", "prim", "tetr", "sexa", "caud", "biau", "mult"]

# Map OMA tags (as they appear in the file) -> your desired label
# NOTE: we do NOT use "burs_hap" as an output label; we map it to "burs" if it ever appears.
TAG_ALIAS = {
    "prim_Ir4": "prim",
    "burs_hap": "burs",
}

# Which raw multi-FASTA to search for each NON-burs species
# (burs is handled separately via per-gene FASTAs)
RAW_PREFIX = {
    "biau": "biau",
    "caud": "caud",
    "mult": "mult",
    "prim": "prim",
    "sexa": "sexa",
    "tetr": "tetr",
}


# ============================================================
# FASTA utilities
# ============================================================

def load_fasta_as_dict(fa_path: Path) -> Dict[str, str]:
    """Read multi-FASTA into dict: {seq_id: sequence}; seq_id is first token after '>'."""
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_seq: List[str] = []

    with fa_path.open("r") as fh:
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


def read_single_fasta_file(fa_path: Path) -> Optional[str]:
    """
    Read a single-record FASTA file and return the sequence string.
    Returns None if file missing or no sequence.
    """
    if not fa_path.exists():
        return None

    seq: List[str] = []
    with fa_path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq.append(line)

    s = "".join(seq)
    return s if s else None


def write_fasta(out_path: Path, records: List[Tuple[str, str]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


# ============================================================
# OMA parsing
# ============================================================

def parse_oma_row(line: str):
    """
    Input is tab-delimited:
      OMAxxxxx <tab> burs:ID <tab> prim_Ir4:ID <tab> tetr:ID ...
    Returns (oma_id, tag2id) where tag2id keys are normalized to KEEP_SPECIES.
    """
    line = line.strip()
    if not line or line.startswith("#"):
        return None, None

    parts = line.split("\t")
    oma_id = parts[0].strip()
    tag2id: Dict[str, str] = {}

    for tok in parts[1:]:
        tok = tok.strip()
        if ":" not in tok:
            continue

        tag, gid = tok.split(":", 1)
        tag = TAG_ALIAS.get(tag.strip(), tag.strip())
        gid = gid.strip()

        if tag in KEEP_SPECIES:
            tag2id.setdefault(tag, gid)  # keep first if duplicates

    return oma_id, tag2id


# ============================================================
# Bursaria per-gene lookup helpers
# ============================================================

def burs_p_to_t(gene_id_p: str) -> str:
    """Convert PBUR.Dd1.1.P####### -> PBUR.Dd1.1.T####### for file lookup."""
    return gene_id_p.replace(".P", ".T", 1)


def fetch_burs_sequence(kind: str, gene_id_p: str) -> Optional[str]:
    """
    kind: 'prot' or 'nucl'
    gene_id_p is from OMA (P-form). Files are named with T-form.
    """
    gene_id_t = burs_p_to_t(gene_id_p)
    fa_dir = BURS_GENE_PROT_DIR if kind == "prot" else BURS_GENE_NUCL_DIR
    fa_path = fa_dir / f"{gene_id_t}.fa"
    return read_single_fasta_file(fa_path)


# ============================================================
# MAIN
# ============================================================

def main():
    # Load non-burs raw multi-FASTAs once
    prot_db: Dict[str, Dict[str, str]] = {}
    nucl_db: Dict[str, Dict[str, str]] = {}

    for sp, prefix in RAW_PREFIX.items():
        prot_fa = RAW_BASE / "protein_seq" / f"{prefix}.fa"
        nucl_fa = RAW_BASE / "nucleotide_seq" / f"{prefix}.fa"

        if not prot_fa.exists():
            raise RuntimeError(f"Missing protein FASTA: {prot_fa}")
        if not nucl_fa.exists():
            raise RuntimeError(f"Missing nucleotide FASTA: {nucl_fa}")

        prot_db[sp] = load_fasta_as_dict(prot_fa)
        nucl_db[sp] = load_fasta_as_dict(nucl_fa)

    # Ensure burs per-gene directories exist
    if not BURS_GENE_PROT_DIR.exists():
        raise RuntimeError(f"Missing burs protein per-gene dir: {BURS_GENE_PROT_DIR}")
    if not BURS_GENE_NUCL_DIR.exists():
        raise RuntimeError(f"Missing burs nucleotide per-gene dir: {BURS_GENE_NUCL_DIR}")

    OUT_PROT.mkdir(parents=True, exist_ok=True)
    OUT_NUCL.mkdir(parents=True, exist_ok=True)
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    wrote = 0
    skipped_missing_tags = 0
    skipped_missing_seqs = 0

    missing_lines: List[str] = ["OMA_ID\tspecies\tgene_id\tmissing_type"]

    with OMA_FILE.open("r") as fh:
        for line in fh:
            oma_id, tag2id = parse_oma_row(line)
            if oma_id is None:
                continue

            total_rows += 1

            # missing species tags from the OMA row (after alias normalization)
            missing_tags = [sp for sp in KEEP_SPECIES if sp not in tag2id]
            for sp in missing_tags:
                missing_lines.append(f"{oma_id}\t{sp}\tNA\tmissing_tag")

            if REQUIRE_ALL_7 and missing_tags:
                skipped_missing_tags += 1
                continue

            prot_records: List[Tuple[str, str]] = []
            nucl_records: List[Tuple[str, str]] = []
            missing_seq_any = False

            for sp in KEEP_SPECIES:
                if sp not in tag2id:
                    continue  # allowed if REQUIRE_ALL_7=False

                gene_id = tag2id[sp]

                # Fetch sequences
                if sp == "burs":
                    pseq = fetch_burs_sequence("prot", gene_id)
                    nseq = fetch_burs_sequence("nucl", gene_id)
                else:
                    pseq = prot_db[sp].get(gene_id)
                    nseq = nucl_db[sp].get(gene_id)

                # Report missingness
                if pseq is None:
                    missing_seq_any = True
                    missing_lines.append(f"{oma_id}\t{sp}\t{gene_id}\tmissing_protein_seq")
                if nseq is None:
                    missing_seq_any = True
                    missing_lines.append(f"{oma_id}\t{sp}\t{gene_id}\tmissing_nucleotide_seq")

                # Add records if present
                if pseq is not None:
                    prot_records.append((f"{sp}|{gene_id}", pseq))
                if nseq is not None:
                    nucl_records.append((f"{sp}|{gene_id}", nseq))

            if REQUIRE_ALL_7 and missing_seq_any:
                skipped_missing_seqs += 1
                continue

            # Write output FASTAs
            write_fasta(OUT_PROT / f"{oma_id}.fa", prot_records)
            write_fasta(OUT_NUCL / f"{oma_id}.fa", nucl_records)
            wrote += 1

    # write report
    with REPORT_TSV.open("w") as rep:
        rep.write("\n".join(missing_lines) + "\n")

    print("Total OMA rows scanned:", total_rows)
    print("FASTA written:", wrote)
    print("Skipped (missing required tags):", skipped_missing_tags)
    print("Skipped (missing sequences in raw FASTAs):", skipped_missing_seqs)
    print("Outputs:")
    print("  Protein   ->", OUT_PROT)
    print("  Nucleotide->", OUT_NUCL)
    print("Missing report ->", REPORT_TSV)


if __name__ == "__main__":
    main()

