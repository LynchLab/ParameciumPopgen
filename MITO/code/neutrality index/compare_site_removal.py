import os, glob, csv



# Root with your per-gene folders containing:

#   - <gene>_sexa_pnps_vs_dnds.fasta

#   - atp_filtered_sites.tsv  (sexa_vs_node_*_REMOVED_codons.tsv copied here)

ROOT = "/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/additonal_species_identification/pnps_sites_0.6/biaurelia/compare_seq_dnds_pnps"



PNPS_HDR_PREFIX = "biaurelia_V1_4"

DNDS_HDR_EXACT  = "biau"



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

                    seq = "".join(buf).replace(" ", "").upper().replace("?", "N")

                    recs.append((h, seq))

                h = line[1:].strip()

                buf = []

            else:

                buf.append(line.strip())

        if h is not None:

            seq = "".join(buf).replace(" ", "").upper().replace("?", "N")

            recs.append((h, seq))

    return recs



def get_seq(recs, exact=None, prefix=None):

    for h,s in recs:

        if exact is not None and h == exact:

            return h,s

        if prefix is not None and h.startswith(prefix):

            return h,s

    return None, None



def to_codons(seq, label):

    L = len(seq)

    if L % 3 != 0:

        raise SystemExit(f"[{label}] length {L} not divisible by 3")

    return [seq[i:i+3] for i in range(0, L, 3)]



def read_removed_indices(tsv_path):

    rm = set()

    with open(tsv_path) as f:

        for ln in f:

            ln = ln.strip()

            if not ln or ln.lower().startswith("codon_index"): 

                continue

            p = ln.split()

            if p and p[0].isdigit():

                rm.add(int(p[0]))

    if not rm:

        raise ValueError(f"No indices parsed from {tsv_path}")

    return rm



def nw_align_tokens(a, b, match=1, mismatch=-1, gap=-1):

    m, n = len(a), len(b)

    score = [[0]*(n+1) for _ in range(m+1)]

    back  = [[None]*(n+1) for _ in range(m+1)]

    for i in range(1,m+1):

        score[i][0] = i*gap; back[i][0] = 'U'

    for j in range(1,n+1):

        score[0][j] = j*gap; back[0][j] = 'L'

    for i in range(1,m+1):

        for j in range(1,n+1):

            s = match if a[i-1]==b[j-1] else mismatch

            diag = score[i-1][j-1] + s

            up   = score[i-1][j]   + gap

            left = score[i][j-1]   + gap

            best, bt = diag, 'D'

            if up > best:   best, bt = up, 'U'

            if left > best: best, bt = left, 'L'

            score[i][j], back[i][j] = best, bt

    A, B = [], []

    i, j = m, n

    while i>0 or j>0:

        step = back[i][j]

        if step=='D':

            A.append(a[i-1]); B.append(b[j-1]); i-=1; j-=1

        elif step=='U':

            A.append(a[i-1]); B.append(None);    i-=1

        else:

            A.append(None);    B.append(b[j-1]); j-=1

    A.reverse(); B.reverse()

    return A, B



ok, skipped = 0, []

for gene in sorted(os.listdir(ROOT)):

    gdir = os.path.join(ROOT, gene)

    if not os.path.isdir(gdir): 

        continue



    # Inputs

    fasta_glob = glob.glob(os.path.join(gdir, "*biau_pnps_vs_dnds.fasta"))

    rm_path    = os.path.join(gdir, "atp_filtered_sites.tsv")

    if not fasta_glob:

        skipped.append((gene, "no *_biau_pnps_vs_dnds.fasta")); continue

    if not os.path.isfile(rm_path):

        skipped.append((gene, "missing atp_filtered_sites.tsv")); continue



    fasta_path = fasta_glob[0]

    try:

        recs = read_fasta(fasta_path)

        _, pnps_seq = get_seq(recs, prefix=PNPS_HDR_PREFIX)

        _, dnds_seq = get_seq(recs, exact=DNDS_HDR_EXACT)

        if pnps_seq is None or dnds_seq is None:

            skipped.append((gene, "missing required headers (biaurelia_V1_4 / biau)")); 

            continue



        pnps_cod = to_codons(pnps_seq, f"{gene}:biaurelia_V1_4")

        dnds_cod = to_codons(dnds_seq, f"{gene}:biau")

        removed  = read_removed_indices(rm_path)



        # If identical sequences, mapping is identity; else do codon NW liftover

        if pnps_seq == dnds_seq:

            pnps_remove = {i for i in removed if 1 <= i <= len(pnps_cod)}

            mapping = {i:i for i in range(1, len(pnps_cod)+1)}

            unmapped = sorted(i for i in removed if i not in mapping)

        else:

            A, B = nw_align_tokens(dnds_cod, pnps_cod)

            mapping = {}

            i_d = i_p = 0

            for a_tok, b_tok in zip(A, B):

                if a_tok is not None: i_d += 1

                if b_tok is not None: i_p += 1

                if a_tok is not None and b_tok is not None:

                    mapping[i_d] = i_p

            pnps_remove = set()

            unmapped = []

            for i in sorted(removed):

                j = mapping.get(i)

                if j is None or not (1 <= j <= len(pnps_cod)):

                    unmapped.append(i)

                else:

                    pnps_remove.add(j)



        # Write output table with PN/PS positions

        out_path = os.path.join(gdir, "biaurelia_V1_4_codon_removal_positions.tsv")

        with open(out_path, "w", newline="") as out:

            w = csv.writer(out, delimiter="\t")

            w.writerow(["pnps_index", "pnps_codon", "dnds_index", "dnds_codon", "note"])

            for i in sorted(removed):

                j = mapping.get(i)

                if j is None or not (1 <= j <= len(pnps_cod)):

                    dcod = dnds_cod[i-1] if 1 <= i <= len(dnds_cod) else ""

                    w.writerow(["", "", i, dcod, "unmapped"])

                else:

                    w.writerow([j, pnps_cod[j-1], i, dnds_cod[i-1], "mapped"])



        print(f"[ok] {gene}: wrote {out_path}  (mapped {len(pnps_remove)} / {len(removed)}; unmapped {len(unmapped)})")

        ok += 1



    except Exception as e:

        skipped.append((gene, f"error: {e}"))



print(f"\nDone. Generated {ok} files.")

if skipped:

    print("\nSkipped:")

    for g, why in skipped:

        print(f"  - {g}: {why}")
