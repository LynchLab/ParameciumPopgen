import os, glob



# --- Paths (as per your layout) ---

BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification"

PNPS_DIR = f"{BASE}/pnps_sites_0.6/biaurelia/codon_alignment"

DNDS_DIR = f"{BASE}/dn_ds/ancestral_dnds"

OUT_ROOT = f"{BASE}/pnps_sites_0.6/biaurelia/compare_seq_dnds_pnps"



PNPS_SEXA_PREFIX = "biaurelia_V1_4"

DNDS_SEXA_HEADER = "biau"



os.makedirs(OUT_ROOT, exist_ok=True)



def read_fasta(path):

    recs = []

    h, buf = None, []

    with open(path) as f:

        for line in f:

            line=line.rstrip("\n")

            if not line:

                continue

            if line.startswith(">"):

                if h is not None:

                    seq = "".join(buf).replace(" ", "")

                    recs.append((h, seq))

                h = line[1:].strip()

                buf = []

            else:

                buf.append(line.strip())

        if h is not None:

            seq = "".join(buf).replace(" ", "")

            recs.append((h, seq))

    return recs



def write_fasta(recs, path, width=70):

    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w") as out:

        for h,s in recs:

            out.write(f">{h}\n")

            for i in range(0, len(s), width):

                out.write(s[i:i+width] + "\n")



def get_by_header_exact_or_prefix(recs, exact=None, prefix=None):

    for h,s in recs:

        if exact is not None and h == exact:

            return h,s

        if prefix is not None and h.startswith(prefix):

            return h,s

    return None, None



success, missing = 0, []



# Loop over every PN/PS codon-aligned file

for fp in sorted(glob.glob(os.path.join(PNPS_DIR, "*_codon_aligned.fna"))):

    gene = os.path.basename(fp).replace("_codon_aligned.fna", "")

    dnds_fa = os.path.join(DNDS_DIR, gene, "all_sequences_cleaned_with_gaps.fasta")



    if not os.path.isfile(dnds_fa):

        missing.append((gene, "missing dN/dS FASTA"))

        continue



    # Read PN/PS gene FASTA and pull 11_II_7__sexaurelia

    pnps_recs = read_fasta(fp)

    h_p, s_p = get_by_header_exact_or_prefix(pnps_recs, prefix=PNPS_SEXA_PREFIX)

    if s_p is None:

        missing.append((gene, f"no PN/PS header starting with '{PNPS_SEXA_PREFIX}'"))

        continue



    # Read dN/dS gene FASTA and pull sexa

    dnds_recs = read_fasta(dnds_fa)

    h_d, s_d = get_by_header_exact_or_prefix(dnds_recs, exact=DNDS_SEXA_HEADER)

    if s_d is None:

        missing.append((gene, "no dN/dS header 'biau'"))

        continue



    out_dir = os.path.join(OUT_ROOT, gene)

    out_fa  = os.path.join(out_dir, f"{gene}_biau_pnps_vs_dnds.fasta")

    write_fasta([(h_p, s_p), (h_d, s_d)], out_fa)

    print(f"[ok] {gene} -> {out_fa}")

    success += 1



print(f"\nDone. Wrote {success} gene files into: {OUT_ROOT}")

if missing:

    print("\nSome genes skipped:")

    for g, why in missing:

        print(f"  - {g}: {why}")
