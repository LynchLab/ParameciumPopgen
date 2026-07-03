#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse


def fasta_headers(fa_path: Path) -> list[str]:
    """Return FASTA header names without leading '>' (first token only)."""
    headers: list[str] = []
    with fa_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def check_folder(base_dir: Path, show_header_diffs: bool = True) -> int:
    """Return number of problems found."""
    prot_dir = base_dir / "protein_seq"
    nuc_dir = base_dir / "nucleotide_seq"

    problems = 0

    if not prot_dir.is_dir():
        print(f"❌ Missing protein directory: {prot_dir}")
        return 1
    if not nuc_dir.is_dir():
        print(f"❌ Missing nucleotide directory: {nuc_dir}")
        return 1

    prot_files = sorted(prot_dir.glob("*.fa"))
    nuc_files = sorted(nuc_dir.glob("*.fa"))

    prot_map = {p.name: p for p in prot_files}
    nuc_map = {p.name: p for p in nuc_files}

    all_names = sorted(set(prot_map) | set(nuc_map))

    print(f"\n=== Checking: {base_dir.resolve()} ===")
    print(f"Protein files: {len(prot_files)} | Nucleotide files: {len(nuc_files)}")

    for name in all_names:
        prot_path = prot_map.get(name)
        nuc_path = nuc_map.get(name)

        if prot_path is None:
            problems += 1
            print(f"❌ Missing protein file for {name}")
            continue
        if nuc_path is None:
            problems += 1
            print(f"❌ Missing nucleotide file for {name}")
            continue

        prot_hdrs = fasta_headers(prot_path)
        nuc_hdrs = fasta_headers(nuc_path)

        prot_n = len(prot_hdrs)
        nuc_n = len(nuc_hdrs)

        if prot_n != nuc_n:
            problems += 1
            print(f"❌ COUNT MISMATCH {name}: protein={prot_n} nucleotide={nuc_n}")
        else:
            print(f"✅ OK {name}: n={prot_n}")

        if show_header_diffs:
            prot_set = set(prot_hdrs)
            nuc_set = set(nuc_hdrs)

            only_prot = sorted(prot_set - nuc_set)
            only_nuc = sorted(nuc_set - prot_set)

            if only_prot or only_nuc:
                problems += 1
                print(f"   ❌ HEADER MISMATCH {name}")
                if only_prot:
                    print("     Protein only:", ", ".join(only_prot))
                if only_nuc:
                    print("     Nucleotide only:", ", ".join(only_nuc))

    if problems == 0:
        print("✅ All pairs match in this folder.")
    else:
        print(f"⚠️ Problems found in this folder: {problems}")

    return problems


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Check if protein_seq and nucleotide_seq FASTA record counts/headers match."
    )
    ap.add_argument(
        "--base",
        type=Path,
        default=Path("."),
        help="Folder that contains protein_seq/ and nucleotide_seq/ (default: current dir).",
    )
    ap.add_argument(
        "--no_header_diffs",
        action="store_true",
        help="Only check counts; do not compare header names.",
    )
    args = ap.parse_args()

    base = args.base.resolve()
    show_header_diffs = not args.no_header_diffs

    problems = check_folder(base, show_header_diffs=show_header_diffs)
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()

