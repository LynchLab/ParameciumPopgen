#!/usr/bin/env python3

import os

import re

import sys

import pandas as pd



# === Inputs / Outputs (run from: .../dn_ds/ancestral_dnds_high_prob) ===

MASTER = "pairwise_numof_syn_nonsyn_gene.tsv"   # tab-separated

OUT    = "combined_normalized_internal_nodes_summary.csv"



# ---------- helpers ----------



def normalize_node(x: str) -> str:

    """Ensure node ids are like 'node_28' regardless of input style."""

    s = str(x).strip()

    if not s:

        return s

    if s.lower().startswith("node_"):

        # normalize case

        return f"node_{re.sub(r'^node_', '', s, flags=re.I)}"

    m = re.search(r'(\d+)$', s)

    return f"node_{m.group(1)}" if m else s



def species_key(s: str) -> str:

    """Strip leading digits + underscore from species labels; keep case as-is."""

    return re.sub(r'^\d+_', '', str(s).strip())



def canon_cols(df: pd.DataFrame) -> pd.DataFrame:

    """Lower + replace spaces/hyphens with underscores."""

    df = df.copy()

    df.columns = [c.strip().lower().replace(" ", "_").replace("-", "_") for c in df.columns]

    return df



def load_all_site_summaries() -> pd.DataFrame:

    """

    Walk subdirectories and collect per-gene site summaries:

    prefer anc_dn_ds_cal_3/summary_sites_by_pair.tsv; fallback to summary.csv.

    Returns a normalized dataframe with columns:

      gene, species, species_filterkey, node, synonymous_sites, non_synonymous_sites

    """

    rows = []

    for d in sorted(x for x in os.listdir(".") if os.path.isdir(x) and not x.startswith(".")):

        cand1 = os.path.join(d, "anc_dn_ds_cal_3", "summary_sites_by_pair.tsv")

        cand2 = os.path.join(d, "anc_dn_ds_cal_3", "summary.csv")

        path = cand1 if os.path.isfile(cand1) else (cand2 if os.path.isfile(cand2) else None)

        if not path:

            continue

        try:

            df = pd.read_csv(path, sep="\t" if path.endswith(".tsv") else ",")

        except Exception as e:

            print(f"⚠️  Skip unreadable summary for gene '{d}': {path} ({e})", file=sys.stderr)

            continue

        df = canon_cols(df)



        # Required/derived columns

        # species_filterkey preferred; if missing, derive from 'species'

        if "species_filterkey" not in df.columns:

            if "species" in df.columns:

                df["species_filterkey"] = df["species"].apply(species_key)

            else:

                print(f"⚠️  {d}: summary lacks 'Species'/'Species_FilterKey' — skipping.", file=sys.stderr)

                continue



        # normalize node col name and values

        node_col = None

        for cand in ("node", "node_resolvedid"):

            if cand in df.columns:

                node_col = cand

                break

        if not node_col:

            print(f"⚠️  {d}: summary lacks a 'Node' column — skipping.", file=sys.stderr)

            continue

        df["node"] = df[node_col].apply(normalize_node)



        # site columns: handle both 'non_synonymous' and 'nonsynonymous'

        syn_col = None

        for cand in ("synonymous_sites", "synonymous_site", "syn_sites"):

            if cand in df.columns:

                syn_col = cand

                break

        nonsyn_col = None

        for cand in ("non_synonymous_sites", "nonsynonymous_sites", "non_syn_sites"):

            if cand in df.columns:

                nonsyn_col = cand

                break

        if syn_col is None or nonsyn_col is None:

            print(f"⚠️  {d}: summary missing synonymous/non-synonymous site columns — skipping.", file=sys.stderr)

            continue



        df_out = pd.DataFrame({

            "gene": d,

            "species": df.get("species", pd.Series([""]*len(df))),

            "species_filterkey": df["species_filterkey"],

            "node": df["node"].apply(normalize_node),

            "synonymous_sites": df[syn_col].astype(float),

            "non_synonymous_sites": df[nonsyn_col].astype(float),

        })

        rows.append(df_out)



    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(

        columns=["gene","species","species_filterkey","node","synonymous_sites","non_synonymous_sites"]

    )



# ---------- main ----------



def main():

    if not os.path.isfile(MASTER):

        print(f"ERROR: master file '{MASTER}' not found", file=sys.stderr)

        sys.exit(1)



    try:

        dfm = pd.read_csv(MASTER, sep="\t")

    except Exception as e:

        print(f"ERROR: failed to read '{MASTER}': {e}", file=sys.stderr)

        sys.exit(1)

    dfm = canon_cols(dfm)



    # Ensure required columns exist

    required = ["gene", "node", "dn", "ds", "species"]

    missing = [c for c in required if c not in dfm.columns]

    if missing:

        print(f"ERROR: master is missing columns: {missing}", file=sys.stderr)

        sys.exit(1)



    # Normalize keys

    dfm["node_norm"] = dfm["node"].apply(normalize_node)

    dfm["species_key"] = dfm["species"].apply(species_key)



    # Load all per-gene site summaries

    sites = load_all_site_summaries()

    if sites.empty:

        print("ERROR: no per-gene site summaries found (look for */anc_dn_ds_cal_3/summary_sites_by_pair.tsv).", file=sys.stderr)

        sys.exit(1)



    # Merge master with sites on (gene, species_key, node)

    merged = dfm.merge(

        sites,

        left_on=["gene", "species_key", "node_norm"],

        right_on=["gene", "species_filterkey", "node"],

        how="left",

        suffixes=("", "_sites")

    )



    # Report missing matches

    missing_mask = merged["synonymous_sites"].isna() | merged["non_synonymous_sites"].isna()

    if missing_mask.any():

        missed = merged.loc[missing_mask, ["gene", "species", "node"]]

        uniq = missed.drop_duplicates()

        for _, r in uniq.iterrows():

            print(f"⚠️  No site summary match for gene={r['gene']} species={r['species']} node={normalize_node(r['node'])}", file=sys.stderr)



    # Compute normalized metrics

    def safe_float(x):

        try:

            return float(x)

        except Exception:

            return float("nan")



    merged["dn"] = merged["dn"].apply(safe_float)

    merged["ds"] = merged["ds"].apply(safe_float)



    syn = merged["synonymous_sites"]

    nonsyn = merged["non_synonymous_sites"]



    # Normalized counts

    merged["Normalized_dN"] = merged["dn"] / nonsyn

    merged["Normalized_dS"] = merged["ds"] / syn



    # Normalized ratio

    merged["Normalized_dN/dS"] = merged["Normalized_dN"] / merged["Normalized_dS"]



    # Pretty rounding and Inf handling

    def fmt_ratio(x):

        if pd.isna(x) or pd.isna(float(x)):

            return ""

        try:

            if x == float("inf") or x == float("-inf"):

                return "Inf"

        except Exception:

            pass

        return round(float(x), 6)



    out_cols = list(dfm.columns)  # keep all original master cols

    # ensure we keep the normalized, and site counts

    extra_cols = [

        "synonymous_sites", "non_synonymous_sites",

        "Normalized_dN", "Normalized_dS", "Normalized_dN/dS"

    ]

    for c in extra_cols:

        if c not in out_cols:

            out_cols.append(c)



    out = merged[out_cols].copy()

    # round numeric outputs

    for c in ["synonymous_sites", "non_synonymous_sites", "Normalized_dN", "Normalized_dS"]:

        if c in out.columns:

            out[c] = out[c].apply(lambda v: (round(v, 6) if pd.notna(v) else ""))

    if "Normalized_dN/dS" in out.columns:

        out["Normalized_dN/dS"] = out["Normalized_dN/dS"].apply(fmt_ratio)



    out.to_csv(OUT, index=False)

    print(f"✅ All results saved to: {OUT}")



if __name__ == "__main__":

    main()


