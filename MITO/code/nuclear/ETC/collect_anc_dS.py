#!/usr/bin/env python3

from pathlib import Path
import shutil

ETC_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/"
    "nuclear_etc_ribo/dnds/etc_ribo/ancestral_dnds_icode5/ETC"
)

OUT_DIR = Path(
    "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/"
    "nuclear_etc_ribo/final_etc_ribo/dnds/etc/ancestral_dS_icode5"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

copied = 0
missing = []

for oma_dir in sorted(ETC_DIR.glob("OMA*")):

    if not oma_dir.is_dir():
        continue

    oma = oma_dir.name

    # Correct input filename
    src_file = oma_dir / "dS.csv"

    # Output filename
    out_file = OUT_DIR / f"{oma}_dS.csv"

    if src_file.is_file():
        shutil.copy2(src_file, out_file)
        copied += 1
        print(f"Copied: {src_file} -> {out_file}")
    else:
        missing.append(oma)
        print(f"Missing: {src_file}")

print("\nSummary")
print("-------")
print(f"ETC folder checked: {ETC_DIR}")
print(f"Output folder: {OUT_DIR}")
print(f"Files copied: {copied}")
print(f"Missing dS.csv files: {len(missing)}")

if missing:
    print("\nFolders missing dS.csv:")
    for oma in missing:
        print(oma)
