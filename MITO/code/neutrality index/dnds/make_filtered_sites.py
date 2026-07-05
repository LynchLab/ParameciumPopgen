#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import math
import shutil
from datetime import datetime

# Run from ancestral_dnds_high_prob/
BASE = Path(".")
THRESHOLD = 0.60

# Input files made by get_site_node_tables.py
INNAME = "closest_node_codonprob_by_site.csv"

# Output folder used by pairwise_dnds_node_species.py and get_syn_nonsyn_sites_3.py
OUTDIR = BASE / "filtered_sites"

# If True, missing/blank probability cells are also excluded
EXCLUDE_MISSING = True


def species_key_from_header(colname: str) -> str:
    """
    Header looks like:
        pent-23 Node
        1_pent-23 Node
        mult-31 Node

    Return:
        pent
        mult
    """
    x = colname.strip()
    x = re.sub(r"\s+Node$", "", x)

    m = re.match(r"(.+)-(\d+)$", x)
    if m:
        sp = m.group(1)
    else:
        sp = x

    # remove leading number + underscore, e.g. 1_pent -> pent
    sp = re.sub(r"^\d+_", "", sp)
    return sp


def parse_prob(cell: str):
    """
    Cell looks like:
        GTA 0.876

    Return:
        0.876
    """
    if cell is None:
        return None

    cell = str(cell).strip()
    if not cell:
        return None

    parts = cell.split()
    if len(parts) < 2:
        return None

    try:
        return float(parts[-1])
    except ValueError:
        return None


def backup_existing_filtered_sites():
    if OUTDIR.exists():
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = BASE / f"filtered_sites_backup_before_prob0.60_{stamp}"
        shutil.copytree(OUTDIR, backup)
        print(f"Backed up old filtered_sites/ to {backup}")

        shutil.rmtree(OUTDIR)


def main():
    backup_existing_filtered_sites()
    OUTDIR.mkdir(exist_ok=True)

    summary_rows = []

    gene_dirs = sorted(
        p for p in BASE.iterdir()
        if p.is_dir()
        and not p.name.startswith(".")
        and p.name not in {"filtered_sites", "anc_site_probs"}
        and not p.name.startswith("filtered_sites_backup")
    )

    wrote_any = False

    for gdir in gene_dirs:
        gene = gdir.name
        prob_csv = gdir / INNAME

        if not prob_csv.exists():
            print(f"Skip {gene}: missing {prob_csv}")
            continue

        with prob_csv.open(newline="") as f:
            rdr = csv.DictReader(f)
            if not rdr.fieldnames or "site" not in rdr.fieldnames:
                print(f"Skip {gene}: no site column in {prob_csv}")
                continue

            species_cols = [c for c in rdr.fieldnames if c != "site"]
            lowprob_sites = {
                species_key_from_header(c): set()
                for c in species_cols
            }

            for row in rdr:
                try:
                    site = int(row["site"])
                except Exception:
                    continue

                for col in species_cols:
                    sp = species_key_from_header(col)
                    prob = parse_prob(row.get(col, ""))

                    if prob is None:
                        if EXCLUDE_MISSING:
                            lowprob_sites[sp].add(site)
                        continue

                    # This means: remove low-confidence ancestral reconstructions
                    if prob < THRESHOLD:
                        lowprob_sites[sp].add(site)

        out_txt = OUTDIR / f"{gene}.txt"

        with out_txt.open("w") as out:
            for sp in sorted(lowprob_sites):
                sites = sorted(lowprob_sites[sp])

                if not sites:
                    continue

                site_string = "-".join(str(x) for x in sites)
                out.write(f"{sp} {site_string}\n")

                summary_rows.append({
                    "gene": gene,
                    "species": sp,
                    "num_filtered_sites": len(sites),
                    "filter_file": str(out_txt)
                })

        wrote_any = True
        print(f"Wrote {out_txt}")

    summary_csv = BASE / "filtered_sites_0.60_summary.csv"
    with summary_csv.open("w", newline="") as f:
        fieldnames = ["gene", "species", "num_filtered_sites", "filter_file"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(summary_rows)

    print(f"\nDone. Summary written to {summary_csv}")

    if not wrote_any:
        print("No filtered_sites files were written. Check that each gene folder has closest_node_codonprob_by_site.csv.")


if __name__ == "__main__":
    main()