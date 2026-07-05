#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Set, Iterable

# ------------------------------------------------------------
# INPUTS
# ------------------------------------------------------------
OMA_FILE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/final_etc_ribo/pnps/etc/OrthologousGroups_ETC_only.txt"
)

# Directory where your pnps gene tables live (the ones you showed with head *.csv)
PNPS_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/etc/gene_pnps/gene_pnps"
)

# ------------------------------------------------------------
# OUTPUT (YOU REQUESTED THIS EXACT LOCATION)
# ------------------------------------------------------------
OUT_BASE = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/"
    "Ravina/nuclear_etc_ribo/final_etc_ribo/pnps/etc"
)
OUT_DIR = OUT_BASE / "pnps_only_OMA_ETC_all_species_fixed"
IDS_DIR = OUT_DIR / "ids_used"

OUT_DIR.mkdir(parents=True, exist_ok=True)
IDS_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Which pnps CSV to filter + which OMA tag to use as the ID source
#
# Includes:
# - the 7 “real” species: biaurelia, bursaria, caudatum, multimicronucleatum,
#   primaurelia, sexaurelia, tetraurelia
# - plus your proxy ones: tribursaria, triaurelia, fokini, quindecaurelia
# ------------------------------------------------------------
CSV_TO_OMA_TAG: Dict[str, str] = {
    # --- 7 core species ---
    "biaurelia.csv":          "biau",
    "bursaria.csv":           "burs_hap",   # IMPORTANT: bursaria uses T => choose burs_hap (T)
    "caudatum.csv":           "caud",
    "multimicronucleatum.csv":"mult",
    "primaurelia.csv":        "prim_Ir4",
    "sexaurelia.csv":         "sexa",
    "tetraurelia.csv":        "tetr",

    # --- proxy species you asked earlier ---
    "tribursaria.csv":        "burs_hap",   # tribursaria uses PBUR...T
    "triaurelia.csv":         "prim_Ir4",   # triaurelia uses PPRIMT..., OMA gives PPRIMP...
    "fokini.csv":             "mult",       # fokini uses PMMNT..., OMA gives PMMNP...
    "quindecaurelia.csv":     "biau",       # quindecaurelia uses PBIA...T, OMA gives PBIA...P
}

# robust extraction even if spacing is messy
PAIR_RE = re.compile(r"([A-Za-z0-9_]+):([^\s]+)")

def extract_ids_for_tag(oma_path: Path, tag: str) -> Set[str]:
    ids: Set[str] = set()
    with oma_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            for t, gene_id in PAIR_RE.findall(line):
                if t == tag:
                    ids.add(gene_id)
    return ids

def add_dotP_to_dotT(ids: Set[str]) -> Set[str]:
    """Add converted IDs where '.P' becomes '.T' (keep original too)."""
    out = set(ids)
    for x in ids:
        if ".P" in x and ".T" not in x:
            out.add(x.replace(".P", ".T", 1))
    return out

def add_prefix_swap(ids: Set[str], old: str, new: str) -> Set[str]:
    """Add converted IDs where prefix old->new (keep original too)."""
    out = set(ids)
    for x in ids:
        if x.startswith(old):
            out.add(x.replace(old, new, 1))
    return out

def convert_ids_to_match_pnps(csv_name: str, oma_tag: str, oma_ids: Set[str]) -> Set[str]:
    """
    Convert OMA IDs to match the gene ID format used in the pnps CSV.
    We ALWAYS keep original OMA IDs and also add converted forms.
    """
    keep = set(oma_ids)

    # bursaria + tribursaria: we deliberately use burs_hap which is already T, no conversion needed.
    if csv_name in {"bursaria.csv", "tribursaria.csv"}:
        return keep

    # Species whose pnps uses ".T" but OMA is ".P"
    if csv_name in {"biaurelia.csv", "caudatum.csv", "sexaurelia.csv", "tetraurelia.csv", "quindecaurelia.csv"}:
        keep = add_dotP_to_dotT(keep)

    # primaurelia/triaurelia: pnps uses PPRIMTxxxxx, OMA prim_Ir4 provides PPRIMPxxxxx
    if csv_name in {"primaurelia.csv", "triaurelia.csv"}:
        keep = add_prefix_swap(keep, "PPRIMP", "PPRIMT")

    # multimicronucleatum/fokini: pnps uses PMMNTxxxxx, OMA mult provides PMMNPxxxxx
    if csv_name in {"multimicronucleatum.csv", "fokini.csv"}:
        keep = add_prefix_swap(keep, "PMMNP", "PMMNT")

    return keep

def filter_pnps_csv(in_csv: Path, out_csv: Path, keep_ids: Set[str]) -> int:
    kept = 0
    with in_csv.open("r", newline="", encoding="utf-8", errors="replace") as fin, \
         out_csv.open("w", newline="", encoding="utf-8") as fout:
        r = csv.reader(fin)
        w = csv.writer(fout)

        header = next(r)
        w.writerow(header)

        for row in r:
            if not row:
                continue
            gene = row[0].strip()
            if gene in keep_ids:
                w.writerow(row)
                kept += 1
    return kept

def write_ids(path: Path, ids: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for x in sorted(set(ids)):
            f.write(x + "\n")

def main() -> None:
    if not OMA_FILE.exists():
        raise FileNotFoundError(f"OMA file not found: {OMA_FILE}")
    if not PNPS_DIR.exists():
        raise FileNotFoundError(f"PNPS dir not found: {PNPS_DIR}")

    for csv_name, oma_tag in CSV_TO_OMA_TAG.items():
        in_csv = PNPS_DIR / csv_name
        if not in_csv.exists():
            print(f"⚠️  Missing pnps file, skipping: {in_csv}")
            continue

        oma_ids = extract_ids_for_tag(OMA_FILE, oma_tag)
        write_ids(IDS_DIR / f"{csv_name}.OMA_tag_{oma_tag}.raw_ids.txt", oma_ids)

        keep_ids = convert_ids_to_match_pnps(csv_name, oma_tag, oma_ids)
        write_ids(IDS_DIR / f"{csv_name}.OMA_tag_{oma_tag}.ids_used_for_filtering.txt", keep_ids)

        out_csv = OUT_DIR / f"{in_csv.stem}.pnps_only_OMA_ETC_using_{oma_tag}.csv"
        kept = filter_pnps_csv(in_csv, out_csv, keep_ids)

        print(f"✅ {csv_name:20s}  tag={oma_tag:8s}  kept={kept}")
        print(f"   out: {out_csv}")

    print(f"\nDone.\nOutput dir: {OUT_DIR}\nIDs saved in: {IDS_DIR}\n")

if __name__ == "__main__":
    main()
