#!/usr/bin/env python3

import re
from pathlib import Path
from Bio import SeqIO

# =========================
# SAME PREFIXES + RULES
# =========================
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

re_dotT = re.compile(r"(" + "|".join(map(re.escape, PREFIXES)) + r")\.T(\d+)")
re_pmmnt  = re.compile(r"\bPMMNT(\d+)\b")     # multimicronucleatum
re_pprimt = re.compile(r"\bPPRIMT(\d+)\b")    # primaurelia

def transcript_to_protein(tid: str) -> str:
    """
    Convert transcript ID -> protein ID
    """
    tid2 = re_dotT.sub(r"\1.P\2", tid)
    tid2 = re_pmmnt.sub(r"PMMNP\1", tid2)
    tid2 = re_pprimt.sub(r"PPRIMP\1", tid2)
    return tid2

def is_transcript_id(s: str) -> bool:
    return (".T" in s) or s.startswith("PMMNT") or s.startswith("PPRIMT")

# =========================
# INPUT DIR
# =========================
BASE = Path("nucleotide_seq")  # run from dnds/etc
if not BASE.exists():
    raise SystemExit(f"❌ Cannot find {BASE.resolve()}")

# Toggle: keep original transcript ID in header description
KEEP_TRANSCRIPT_IN_DESC = True

# =========================
# MAIN
# =========================
n_total = 0
n_changed = 0
n_skipped = 0

for fa in sorted(BASE.rglob("*.fa")):
    n_total += 1

    # transcript id typically equals filename stem and record.id
    old_stem = fa.stem

    # parse fasta (expect single record per file)
    recs = list(SeqIO.parse(str(fa), "fasta"))
    if not recs:
        print(f"⚠️ Empty FASTA: {fa}")
        n_skipped += 1
        continue
    if len(recs) != 1:
        print(f"⚠️ Multiple records in {fa} (skipping)")
        n_skipped += 1
        continue

    rec = recs[0]
    old_id = rec.id  # first token after '>'

    # decide which identifier to convert (prefer header id)
    source_id = old_id if old_id else old_stem

    if not is_transcript_id(source_id):
        # already protein-like
        n_skipped += 1
        continue

    new_id = transcript_to_protein(source_id)

    # update FASTA header
    rec.id = new_id
    rec.name = new_id
    if KEEP_TRANSCRIPT_IN_DESC:
        rec.description = f"{new_id} transcript={source_id}"
    else:
        rec.description = new_id

    # write to temp then replace (safe)
    tmp_path = fa.with_suffix(".fa.tmp")
    SeqIO.write(rec, str(tmp_path), "fasta")

    # rename file
    new_path = fa.with_name(new_id + ".fa")

    # if target exists, avoid overwriting
    if new_path.exists() and new_path.resolve() != fa.resolve():
        print(f"⚠️ Target exists, skipping rename: {new_path}")
        tmp_path.unlink(missing_ok=True)
        n_skipped += 1
        continue

    # replace original content, then rename (if needed)
    tmp_path.replace(fa)
    if new_path.resolve() != fa.resolve():
        fa.rename(new_path)

    n_changed += 1

print(f"\n✅ Done. Total: {n_total} | Updated: {n_changed} | Skipped: {n_skipped}")

