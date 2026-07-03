import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# ============= USER SETTINGS =============
IN_PATH  = "gc_content_per_gene.csv"
REMOVE_rRNA_tRNA = True              # set False to keep rRNA/tRNA
NORMALIZE_BY_GENE_MEDIAN = True     # True -> value / per-gene median; False -> raw %
OUT_RAW  = "aurelia_gc_percent_gray_boxplots_colored_others.pdf"
OUT_NORM = "aurelia_gc_normalized_gray_boxplots_colored_others.pdf"
# ========================================

# ---------- Load ----------
df = pd.read_csv(IN_PATH)

# ---------- Standardize column names ----------
# Try to detect the "sample", "gene_name", and "gc" columns even if named differently
col_map_candidates = {
    "sample":     ["sample", "species", "species_name", "sample_name", "strain", "taxon"],
    "gene_name":  ["gene_name", "gene", "locus", "feature"],
    "gc_value":   ["gc_percent", "gc_content", "gc", "GC%", "gc_pct"]
}

def find_col(df_cols, choices):
    lc = {c.lower(): c for c in df_cols}
    for cand in choices:
        if cand.lower() in lc:
            return lc[cand.lower()]
    return None

sample_col   = find_col(df.columns, col_map_candidates["sample"])
gene_col     = find_col(df.columns, col_map_candidates["gene_name"])
gc_val_col   = find_col(df.columns, col_map_candidates["gc_value"])

missing_bits = []
if sample_col is None: missing_bits.append("sample")
if gene_col   is None: missing_bits.append("gene_name")
if gc_val_col is None: missing_bits.append("gc_percent")
if missing_bits:
    raise ValueError(
        "Missing required columns. Could not find: "
        + ", ".join(missing_bits)
        + f"\nAvailable columns: {list(df.columns)}\n"
        "Expected something like: sample/species, gene/gene_name, gc_percent/gc_content/gc"
    )

# Rename to canonical names used below
df = df.rename(columns={
    sample_col: "sample",
    gene_col: "gene_name",
    gc_val_col: "gc_percent"
})

# ---------- Clean ----------
df["sample"] = df["sample"].astype(str).str.strip()
df["gene_name"] = df["gene_name"].astype(str).str.strip()
df["gc_percent"] = pd.to_numeric(df["gc_percent"], errors="coerce")

# If GC is 0–1, convert to %
if df["gc_percent"].max(skipna=True) is not None and df["gc_percent"].max(skipna=True) <= 1.05:
    df["gc_percent"] = df["gc_percent"] * 100.0

# ---------- Remove tRNA / rRNA (optional) ----------
if REMOVE_rRNA_tRNA:
    mask = ~df["gene_name"].str.lower().str.contains(r"trna|rrna|rna")
    df = df[mask].copy()

# ---------- Optional normalization (value / per-gene median) ----------
if NORMALIZE_BY_GENE_MEDIAN:
    df["yval"] = df.groupby("gene_name")["gc_percent"].transform(lambda x: x / x.median())
    y_col = "yval"
    y_label = "GC% / Median (per gene)"
    out_name = OUT_NORM
else:
    y_col = "gc_percent"
    y_label = "GC content (%)"
    out_name = OUT_RAW

# ---------- Define species groups ----------
# Use the first 4 letters of sample as a key (matches your prior abbreviations)
def abbrev(s):
    return str(s).strip().lower()[:4]

df["group_key"] = df["sample"].map(abbrev)

aurelia_keys = {
    "tetr", "tria", "pent", "sexa", "octa", "dode", "tred",
    "biau", "deca", "jenn", "nova", "prim", "psch", "quad", "sonn"
}
others_keys = {"caud", "mult", "foki", "burs", "tbur", "calk"}

subset = df[df["group_key"].isin(aurelia_keys | others_keys)].copy()

aurelia_df = subset[subset["group_key"].isin(aurelia_keys)].copy()
others_df  = subset[subset["group_key"].isin(others_keys)].copy()

# ---------- Gene order ----------
gene_order = sorted(subset["gene_name"].unique())

# ---------- Plot ----------
plt.figure(figsize=(22, 8))

# Aurelia boxplots (gray)
sns.boxplot(
    data=aurelia_df,
    x="gene_name",
    y=y_col,
    order=gene_order,
    color="lightgray",
    showfliers=False,
    zorder=1
)

# Aurelia points (gray)
sns.stripplot(
    data=aurelia_df,
    x="gene_name",
    y=y_col,
    order=gene_order,
    color="gray",
    size=4,
    jitter=0.18,
    alpha=0.7,
    zorder=2
)

# Non-Aurelia colored points
# Build a stable hue order for the 'others' by sample label (not just the 4-letter key)
others_hue_order = (
    others_df.groupby("sample", as_index=False)
    .size()
    .sort_values("sample")["sample"]
    .tolist()
)
palette = sns.color_palette("tab10", n_colors=max(1, len(others_hue_order)))

sns.stripplot(
    data=others_df,
    x="gene_name",
    y=y_col,
    order=gene_order,
    hue="sample",
    hue_order=others_hue_order,
    dodge=True,
    size=5,
    jitter=0.18,
    palette=palette,
    zorder=3
)

# ---------- Styling ----------
title_suffix = " (normalized by per-gene median)" if NORMALIZE_BY_GENE_MEDIAN else " (raw %)"
plt.title("GC Content per Gene: Aurelia (Gray) vs Non-Aurelia (Colored)" + title_suffix,
          fontsize=18, weight="bold")
plt.xlabel("Gene", fontsize=13)
plt.ylabel(y_label, fontsize=13)
plt.xticks(rotation=90, fontsize=11)
plt.yticks(fontsize=11)

# Show legend only if we have others
if not others_df.empty:
    plt.legend(title="Non-Aurelia Species", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
else:
    leg = plt.legend()
    if leg:
        leg.remove()

plt.tight_layout()

# ---------- Save ----------
plt.savefig(out_name)
plt.close()

print(f"✅ Saved: {out_name}")
print("Columns used:", list(df.columns))
