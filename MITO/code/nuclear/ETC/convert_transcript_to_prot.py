#!/usr/bin/env python3

import re
from pathlib import Path
import pandas as pd

# ==================================================
# 1) SAME PREFIXES (exact mirror of your script)
# ==================================================

PREFIXES = [
    "PBIA.V1_4.1",
    "PBUR.Dd1.1",
    "PCAU.43c3d.1",
    "PDEC.223.1",
    "PDODEC.274.1",
    "PJENN.M.1",
    "PNOV.TE.1",
    "POCT.K8.1",
    "PPENT.87.1",
    "PQUADEC.NiA.1",
    "PSEX.AZ8_4.1",
    "PSON.ATCC_30995.1",
    "PTET.51.1",
    "PTRED.209.2",
]

# .T##### → .P#####
re_dotT = re.compile(
    r'(' + '|'.join(map(re.escape, PREFIXES)) + r')\.T(\d+)'
)

# special cases
re_pmmnt  = re.compile(r'\bPMMNT(\d+)\b')    # multimicronucleatum
re_pprimt = re.compile(r'\bPPRIMT(\d+)\b')   # primaurelia


def reverse_id(tid: str) -> str:
    """
    Convert transcript ID → protein ID
    """
    tid = re_dotT.sub(r'\1.P\2', tid)
    tid = re_pmmnt.sub(r'PMMNP\1', tid)
    tid = re_pprimt.sub(r'PPRIMP\1', tid)
    return tid


# ==================================================
# 2) Paths
# ==================================================

IN_DIR  = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc/id")
OUT_DIR = IN_DIR / "prot_id"

OUT_DIR.mkdir(exist_ok=True)

# ==================================================
# 3) Convert all CSV files
# ==================================================

files = sorted(IN_DIR.glob("*.csv"))
if not files:
    raise SystemExit("❌ No .csv files found in input directory")

for infile in files:
    df = pd.read_csv(infile)

    if "gene" not in df.columns:
        print(f"⚠️ Skipping {infile.name} (no gene column)")
        continue

    df["gene"] = df["gene"].astype(str).apply(reverse_id)

    out_fp = OUT_DIR / infile.name
    df.to_csv(out_fp, index=False)

    print(f"✅ {infile.name} → prot_id/{infile.name}")

print("\n🎉 All transcript IDs converted back to protein IDs")

