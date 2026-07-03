#!/usr/bin/env python3

import os

import sys

import argparse

import shutil

import pandas as pd

from Bio import SeqIO





def fail(msg):

    print(f"ERROR: {msg}", file=sys.stderr)

    sys.exit(1)





def read_fasta_dict(p):

    if not os.path.isfile(p):

        fail(f"Missing FASTA: {p}")

    return {r.id: r for r in SeqIO.parse(p, "fasta")}





def pick_blast_file(blast_dir, sample):

    # Prefer with_coords, then filtered, then raw

    cands = [

        os.path.join(blast_dir, f"{sample}_blastp_with_coords.tsv"),

        os.path.join(blast_dir, f"{sample}_blastp_filtered.tsv"),

        os.path.join(blast_dir, f"{sample}_blastp.tsv"),

    ]

    for p in cands:

        if os.path.isfile(p) and os.path.getsize(p) > 0:

            return p

    return None





def load_annot_intervals(blast_path, coords_df):

    """Return DataFrame with columns start, stop (ints) for annotated intervals."""

    if not blast_path:

        return pd.DataFrame(columns=["start", "stop"])



    # try headered first

    try:

        df = pd.read_csv(blast_path, sep="\t", header=0, dtype=str)

        if "qseqid" not in df.columns and "query" in df.columns:

            df = df.rename(columns={"query": "qseqid"})

    except Exception:

        try:

            df = pd.read_csv(blast_path, sep="\t", header=None, dtype=str)

        except pd.errors.EmptyDataError:

            return pd.DataFrame(columns=["start", "stop"])



    # If start/stop present, use directly

    if "start" in df.columns and "stop" in df.columns:

        out = df[["start", "stop"]].dropna().copy()

        out["start"] = pd.to_numeric(out["start"], errors="coerce")

        out["stop"] = pd.to_numeric(out["stop"], errors="coerce")

        out = out.dropna()

        out["start"] = out["start"].astype(int)

        out["stop"] = out["stop"].astype(int)

        return out



    # Else reconstruct via qseqid -> cds/frame -> coords

    if "qseqid" not in df.columns:

        if df.shape[1] == 0:

            return pd.DataFrame(columns=["start", "stop"])

        df = df.rename(columns={df.columns[0]: "qseqid"})



    jj = df[["qseqid"]].copy()

    jj["cds"] = jj["qseqid"].str.extract(r"^cds_(\d+)_", expand=False).astype("Int64")

    jj["frame"] = jj["qseqid"].str.extract(r"_frame_(\d+)$", expand=False).astype("Int64")

    jj = jj.dropna(subset=["cds", "frame"])



    if jj.empty or coords_df.empty or not {"cds", "frame", "start", "stop"}.issubset(coords_df.columns):

        return pd.DataFrame(columns=["start", "stop"])



    merged = jj.merge(coords_df, on=["cds", "frame"], how="left")

    out = merged[["start", "stop"]].dropna().copy()

    out["start"] = pd.to_numeric(out["start"], errors="coerce")

    out["stop"] = pd.to_numeric(out["stop"], errors="coerce")

    out = out.dropna()

    out["start"] = out["start"].astype(int)

    out["stop"] = out["stop"].astype(int)

    return out





def merge_intervals(df):

    """Merge overlapping/adjacent intervals for speed."""

    if df.empty:

        return []

    ivals = [

        (min(int(a), int(b)), max(int(a), int(b)))

        for a, b in df[["start", "stop"]].itertuples(index=False)

    ]

    ivals.sort()

    merged = []

    for s, e in ivals:

        if not merged or s > merged[-1][1] + 1:

            merged.append([s, e])

        else:

            merged[-1][1] = max(merged[-1][1], e)

    return [(s, e) for s, e in merged]





def overlaps(a, b, c, d):

    """Strand-agnostic interval overlap."""

    a, b = (a, b) if a <= b else (b, a)

    c, d = (c, d) if c <= d else (d, c)

    return not (b < c or d < a)





def process_sample(sample_dir, out_root):

    sample_dir = os.path.abspath(sample_dir)

    sample = os.path.basename(sample_dir)

    species_dir = os.path.basename(os.path.dirname(sample_dir))  # e.g. biaurelia_fasta

    species = species_dir.replace("_fasta", "")



    prot_faa = os.path.join(sample_dir, "protein.faa")

    nuc_fna = os.path.join(sample_dir, "cds_from_genomic.fna")

    coords_csv = os.path.join(sample_dir, "cds_coords.csv")

    blast_dir = os.path.join(sample_dir, "blastp")



    if not os.path.isdir(sample_dir):

        fail(f"Sample dir not found: {sample_dir}")

    if not os.path.isfile(prot_faa):

        fail(f"Missing protein.faa: {prot_faa}")

    if not os.path.isfile(coords_csv):

        fail(f"Missing coords: {coords_csv}")

    if not os.path.isdir(blast_dir):

        fail(f"Missing blastp dir: {blast_dir}")



    coords = pd.read_csv(coords_csv)

    for c in ["cds", "frame", "start", "stop"]:

        if c not in coords.columns:

            fail(f"coords missing column: {c}")

    coords["cds"] = pd.to_numeric(coords["cds"], errors="coerce").astype("Int64")

    coords["frame"] = pd.to_numeric(coords["frame"], errors="coerce").astype("Int64")

    coords["start"] = pd.to_numeric(coords["start"], errors="coerce").astype("Int64")

    coords["stop"]  = pd.to_numeric(coords["stop"],  errors="coerce").astype("Int64")

    coords = coords.dropna(subset=["cds", "frame", "start", "stop"]).astype(

        {"cds": int, "frame": int, "start": int, "stop": int}

    )

    coords["orf_id"] = coords.apply(

        lambda r: f"cds_{r.cds}_frame_{r.frame}", axis=1

    )



    # Load protein / nucleotide sequences

    prot = read_fasta_dict(prot_faa)

    nuc = read_fasta_dict(nuc_fna)



    # Choose annotated BLAST file (priority with_coords)

    annot_file = pick_blast_file(blast_dir, sample)

    if annot_file:

        print(f"[INFO] {species}/{sample}: Annotated source: {annot_file}")

    else:

        print(f"[WARN] {species}/{sample}: No annotated file found; treating as zero annotated.")



    # get annotated intervals (for overlap)

    ann_df = load_annot_intervals(annot_file, coords)

    merged = merge_intervals(ann_df)



    # ------------------ annotated set (from BLAST) ------------------

    bwc = os.path.join(blast_dir, f"{sample}_blastp_with_coords.tsv")

    blast_for_ann = bwc if os.path.isfile(bwc) else annot_file



    ann_orf_ids = set()

    df_ann = None

    if blast_for_ann and os.path.isfile(blast_for_ann):

        df_ann = pd.read_csv(blast_for_ann, sep="\t", dtype=str)

        if "qseqid" in df_ann.columns:

            qcol = "qseqid"

        elif "query" in df_ann.columns:

            qcol = "query"

        else:

            qcol = df_ann.columns[0]

            df_ann = df_ann.rename(columns={qcol: "qseqid"})

            qcol = "qseqid"

        ann_orf_ids = set(df_ann[qcol].dropna().unique())

    else:

        print(f"[WARN] {species}/{sample}: No BLAST file for annotated output.")



    # ------------------ unannotated candidates ------------------

    # start from coords that are NOT in BLAST (not annotated)

    unann = coords[~coords["orf_id"].isin(ann_orf_ids)].copy()



    # filter by non-overlap with ANY annotated interval

    if merged:

        keep_rows = []

        for _, r in unann.iterrows():

            s, e = int(r.start), int(r.stop)

            ov = any(overlaps(s, e, ms, me) for (ms, me) in merged)

            if not ov:

                keep_rows.append(r)

        unann = pd.DataFrame(keep_rows) if keep_rows else pd.DataFrame(columns=coords.columns)

    # else: no annotated intervals → all non-BLAST ORFs are intergenic candidates



    # add aa_len and keep only ≥100 aa

    if not unann.empty:

        def aa_len(orf_id):

            rec = prot.get(orf_id)

            return len(rec.seq) if rec is not None else None



        unann["aa_len"] = unann["orf_id"].map(aa_len)

        unann = unann.dropna(subset=["aa_len"])

        unann = unann[unann["aa_len"] >= 100].copy()



    print(f"[INFO] {species}/{sample}: Unannotated non-overlap >=100 aa: {len(unann)}")



    # ------------------ write outputs under intergenic_region ------------------

    base_out = os.path.join(out_root, species, sample)

    ann_out = os.path.join(base_out, "protein_annotated")

    unann_out = os.path.join(base_out, "protein_unannotated")

    os.makedirs(ann_out, exist_ok=True)

    os.makedirs(unann_out, exist_ok=True)



    prot_dict = prot

    nuc_dict = nuc



    # ========== annotated FASTAs ==========

    if df_ann is not None:

        if "qseqid" in df_ann.columns:

            qcol = "qseqid"

        elif "query" in df_ann.columns:

            qcol = "query"

        else:

            qcol = df_ann.columns[0]

            df_ann = df_ann.rename(columns={qcol: "qseqid"})

            qcol = "qseqid"



        if "subject" in df_ann.columns:

            scol = "subject"

        elif "sseqid" in df_ann.columns:

            scol = "sseqid"

        else:

            scol = None



        ann_prot_out = os.path.join(ann_out, "annotated_protein.faa")

        ann_nuc_out = os.path.join(ann_out, "annotated_nuc.fna")



        ann_prot_records = []

        ann_nuc_records = []



        for _, row in df_ann.dropna(subset=[qcol]).iterrows():

            oid = row[qcol]

            if oid not in prot_dict or oid not in nuc_dict:

                continue

            prot_rec = prot_dict[oid]

            nuc_rec = nuc_dict[oid]



            if scol and isinstance(row.get(scol, None), str):

                subj = row[scol]

                gene_name = subj.split("_")[-1] if "_" in subj else subj

            else:

                gene_name = oid



            new_id = gene_name

            prot_rec.id = new_id

            prot_rec.description = f"{new_id} ({oid})"

            nuc_rec.id = new_id

            nuc_rec.description = f"{new_id} ({oid})"



            ann_prot_records.append(prot_rec)

            ann_nuc_records.append(nuc_rec)



        with open(ann_prot_out, "w") as fp:

            SeqIO.write(ann_prot_records, fp, "fasta")

        with open(ann_nuc_out, "w") as fp:

            SeqIO.write(ann_nuc_records, fp, "fasta")



        # copy blastp_with_coords (if exists) or annot_file into protein_annotated

        src_blast = bwc if os.path.isfile(bwc) else annot_file

        if src_blast:

            shutil.copy2(src_blast, os.path.join(ann_out, os.path.basename(src_blast)))

    else:

        print(f"[WARN] {species}/{sample}: Skipping annotated FASTA (no BLAST).")



    # ========== unannotated FASTAs & CSV ==========

    if not unann.empty:

        unann_ids = sorted(unann["orf_id"].dropna().unique())

        unann_prot_out = os.path.join(unann_out, "unannotated_protein.faa")

        unann_nuc_out = os.path.join(unann_out, "unannotated_nuc.fna")

        unann_csv_out = os.path.join(unann_out, "unannotated.csv")



        unann_prot_records = []

        unann_nuc_records = []



        for oid in unann_ids:

            if oid in prot_dict:

                r = prot_dict[oid]

                r.id = oid

                r.description = oid

                unann_prot_records.append(r)

            if oid in nuc_dict:

                r = nuc_dict[oid]

                r.id = oid

                r.description = oid

                unann_nuc_records.append(r)



        with open(unann_prot_out, "w") as fp:

            SeqIO.write(unann_prot_records, fp, "fasta")

        with open(unann_nuc_out, "w") as fp:

            SeqIO.write(unann_nuc_records, fp, "fasta")



        # keep aa_len and coords in CSV

        cols = ["orf_id", "cds", "frame", "start", "stop", "aa_len"]

        keep_cols = [c for c in cols if c in unann.columns]

        unann[keep_cols].to_csv(unann_csv_out, sep="\t", index=False)

    else:

        print(f"[INFO] {species}/{sample}: No unannotated non-overlap >=100 aa to write.")



    print(f"[DONE] {species}/{sample}: organized into intergenic_region")





def main():

    ap = argparse.ArgumentParser(

        description="Organize annotated and unannotated intergenic ORFs for all samples."

    )

    ap.add_argument(

        "--base-dir",

        required=True,

        help="Root directory containing *_fasta species dirs (e.g., additional_orf)",

    )

    ap.add_argument(

        "--out-root",

        required=True,

        help="Output root for intergenic_region (species/sample/protein_annotated|protein_unannotated)",

    )

    args = ap.parse_args()



    base_dir = os.path.abspath(args.base_dir)

    out_root = os.path.abspath(args.out_root)



    if not os.path.isdir(base_dir):

        fail(f"Base dir not found: {base_dir}")



    for species_dir in sorted(os.listdir(base_dir)):

        s_path = os.path.join(base_dir, species_dir)

        if not os.path.isdir(s_path):

            continue

        if not species_dir.endswith("_fasta"):

            continue

        print(f"=== Species: {species_dir} ===")

        for sample_dir in sorted(os.listdir(s_path)):

            if not sample_dir.endswith("_orf"):

                continue

            full_sample = os.path.join(s_path, sample_dir)

            if not os.path.isdir(full_sample):

                continue

            print(f"  -> Sample: {sample_dir}")

            try:

                process_sample(full_sample, out_root)

            except Exception as e:

                print(f"[ERROR] Failed on {species_dir}/{sample_dir}: {e}", file=sys.stderr)





if __name__ == "__main__":

    main()


