import os, glob, shutil



BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification"

DNDS = f"{BASE}/dn_ds/ancestral_dnds_filterallsites"

DEST = f"{BASE}/pnps_sites_0.6/biaurelia/compare_seq_dnds_pnps"



ok, skipped = 0, []



for gene in sorted(os.listdir(DNDS)):

    gdir = os.path.join(DNDS, gene)

    if not os.path.isdir(gdir): 

        continue

    indd = os.path.join(gdir, "internal_node_dnds")

    if not os.path.isdir(indd):

        continue

    # Prefer sexa_ pattern, fallback to sexaurelia_

    paths = glob.glob(os.path.join(indd, "biau_vs_node_*_REMOVED_codons.tsv"))

    if not paths:

        paths = glob.glob(os.path.join(indd, "biaurelia_vs_node_*_REMOVED_codons.tsv"))

    if not paths:

        skipped.append((gene, "no biau*_REMOVED_codons.tsv"))

        continue

    src = sorted(paths)[0]

    outdir = os.path.join(DEST, gene)

    os.makedirs(outdir, exist_ok=True)

    dst = os.path.join(outdir, "atp_filtered_sites.tsv")

    shutil.copyfile(src, dst)

    print(f"[ok] {gene} -> {dst}  (src: {os.path.basename(src)})")

    ok += 1



print(f"\nDone. Copied {ok} files into {DEST}")

if skipped:

    print("\nSkipped:")

    for g, why in skipped:

        print(f"  - {g}: {why}")
