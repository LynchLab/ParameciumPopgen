import os, glob



BASE = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification"

SEXA = f"{BASE}/pnps_sites_0.6/biaurelia"

ALIGN_IN  = f"{SEXA}/codon_alignment"

ALIGN_OUT = f"{SEXA}/codon_alignment_filtered"

COMPARE   = f"{SEXA}/compare_seq_dnds_pnps"



os.makedirs(ALIGN_OUT, exist_ok=True)



def read_fasta(path):

    recs = []

    h, buf = None, []

    with open(path) as f:

        for line in f:

            line = line.rstrip("\n")

            if not line: 

                continue

            if line.startswith(">"):

                if h is not None:

                    recs.append((h, "".join(buf).replace(" ", "").upper()))

                h = line[1:].strip()

                buf = []

            else:

                buf.append(line.strip())

        if h is not None:

            recs.append((h, "".join(buf).replace(" ", "").upper()))

    return recs



def write_fasta(records, path, width=70):

    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, "w") as out:

        for h, s in records:

            out.write(f">{h}\n")

            for i in range(0, len(s), width):

                out.write(s[i:i+width] + "\n")



def to_codons(seq, label):

    if len(seq) % 3 != 0:

        raise ValueError(f"[{label}] length {len(seq)} not divisible by 3")

    return [seq[i:i+3] for i in range(0, len(seq), 3)]



def parse_pnps_indices(tsv_path):

    # Robustly find the 'pnps_index' column (tab/space separated)

    with open(tsv_path) as f:

        header = None

        idx_col = None

        indices = set()

        for ln in f:

            ln = ln.strip()

            if not ln:

                continue

            if header is None:

                header = ln.split()

                # Locate 'pnps_index' column

                try:

                    idx_col = header.index("pnps_index")

                except ValueError:

                    # fallback: first column

                    idx_col = 0

                continue

            parts = ln.split()

            if idx_col < len(parts):

                val = parts[idx_col]

                # accept pure integers only (skip blanks/unmapped)

                if val.isdigit():

                    i = int(val)

                    if i > 0:

                        indices.add(i)

        if not indices:

            raise ValueError(f"No PN/PS indices parsed from {tsv_path}")

        return indices



def filter_alignment_by_indices(records, rm_1based, gene):

    # Sanity: all sequences same length & multiple of 3

    lengths = {len(s) for _, s in records}

    if len(lengths) != 1:

        raise ValueError(f"[{gene}] sequences not equal length: {sorted(lengths)}")

    L = next(iter(lengths))

    if L % 3 != 0:

        raise ValueError(f"[{gene}] alignment length {L} not multiple of 3")

    n_codons = L // 3



    # Bound check

    bad = [i for i in rm_1based if not (1 <= i <= n_codons)]

    if bad:

        # Keep only valid ones; warn via print

        print(f"[warn] {gene}: {len(bad)} indices out of range (1..{n_codons}) ignored: {bad[:10]}")



    rm = {i for i in rm_1based if 1 <= i <= n_codons}



    # Apply

    out = []

    for h, s in records:

        cods = [s[i*3:(i+1)*3] for i in range(n_codons)]

        kept = [c for k, c in enumerate(cods, start=1) if k not in rm]

        out.append((h, "".join(kept)))

    return out, n_codons, len(rm)



ok, sk = 0, []

for gene_dir in sorted(glob.glob(os.path.join(COMPARE, "*"))):

    if not os.path.isdir(gene_dir):

        continue

    gene = os.path.basename(gene_dir)

    tsv = os.path.join(gene_dir, "biaurelia_V1_4_codon_removal_positions.tsv")

    fin = os.path.join(ALIGN_IN, f"{gene}_codon_aligned.fna")

    fout = os.path.join(ALIGN_OUT, f"{gene}_codon_aligned_filtered.fna")



    if not os.path.isfile(tsv):

        sk.append((gene, "missing biaurelia_V1_4_codon_removal_positions.tsv")); continue

    if not os.path.isfile(fin):

        sk.append((gene, "missing input alignment")); continue



    try:

        rm_idx = parse_pnps_indices(tsv)

        recs   = read_fasta(fin)

        filt, total_cod, removed = filter_alignment_by_indices(recs, rm_idx, gene)

        write_fasta(filt, fout)

        kept = total_cod - removed

        print(f"[ok] {gene}: codons total={total_cod}, removed={removed}, kept={kept} -> {fout}")

        ok += 1

    except Exception as e:

        sk.append((gene, f"error: {e}"))



print(f"\nDone. Filtered {ok} genes into: {ALIGN_OUT}")

if sk:

    print("\nSkipped:")

    for g, why in sk:

        print(f"  - {g}: {why}")
