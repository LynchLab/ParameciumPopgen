import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =========================================================
# INPUT FILES
# =========================================================
gene_length_file = "combined_gene_coords_table.csv"
gc_file = "gc_content_per_gene.csv"

# =========================================================
# SETTINGS
# =========================================================
OUT_FIG = "gene_length_and_gc_all_genes_stacked.png"
REMOVE_RRNA_TRNA = True
NORMALIZE_GENE_LENGTH = True   # keep gene length normalized
USE_RAW_GC = True              # raw GC (no normalization)

# Aurelia and outgroup sample codes
aurelia = [
    "tetr", "tria", "pent", "sexa", "octa", "dode", "tred",
    "biau", "deca", "jenn", "nova", "prima", "psch", "quad", "sonn"
]
others = ["caud", "mult", "foki", "burs", "tbur", "calk"]

# Map sample codes to display names
sample_name_map = {
    "caud": "P. caudatum",
    "burs": "P. bursaria",
    "tbur": "P. tribursaria",
    "mult": "P. multimicronucleatum",
    "foki": "P. fokini",
    "calk": "P. calkinsi",
}

# Outgroup colors
outgroup_colors = {
    "P. tribursaria": "#88c96b",
    "P. bursaria": "#cc0000",
    "P. multimicronucleatum": "#9874f4",
    "P. fokini": "#ffe599",
    "P. caudatum": "#e69138",
    "P. calkinsi": "#000000",
}

# =========================================================
# GLOBAL STYLE
# =========================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 24,          # 🔥 increase base font
    "axes.labelsize": 24,     # 🔥 axis labels bigger
    "axes.titlesize": 30,     # 🔥 titles if you add them
    "xtick.labelsize": 24,    # 🔥 x tick labels
    "ytick.labelsize": 24,    # 🔥 y tick labels
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "text.color": "black",
    "axes.edgecolor": "black",
})

# =========================================================
# HELPER
# =========================================================
def abbrev(s):
    return str(s).strip().lower()[:4]

# =========================================================
# LOAD GENE LENGTH DATA
# =========================================================
len_df = pd.read_csv(gene_length_file)

len_df["sample"] = len_df["sample"].astype(str).str.strip()
len_df["gene_name"] = len_df["gene_name"].astype(str).str.strip()

if REMOVE_RRNA_TRNA:
    len_df = len_df[
        ~len_df["gene_name"].str.lower().str.contains(r"trna|rrna|rna", regex=True)
    ].copy()

if NORMALIZE_GENE_LENGTH:
    len_df["plot_length"] = len_df.groupby("gene_name")["length"].transform(
        lambda x: x / x.median()
    )
    len_ylabel = "Normalized gene length"
else:
    len_df["plot_length"] = len_df["length"]
    len_ylabel = "Gene length (bp)"

len_df["species"] = len_df["sample"].map(sample_name_map).fillna(len_df["sample"])
len_aurelia_df = len_df[len_df["sample"].isin(aurelia)].copy()
len_others_df  = len_df[len_df["sample"].isin(others)].copy()

# =========================================================
# LOAD GC DATA
# =========================================================
gc_df = pd.read_csv(gc_file)

col_map_candidates = {
    "sample": ["sample", "species", "species_name", "sample_name", "strain", "taxon"],
    "gene":   ["gene_name", "gene", "locus", "feature"],
    "gc":     ["gc_percent", "gc_content", "gc", "GC%", "gc_pct"]
}

def find_col(df_cols, choices):
    lc = {c.lower(): c for c in df_cols}
    for cand in choices:
        if cand.lower() in lc:
            return lc[cand.lower()]
    return None

sample_col = find_col(gc_df.columns, col_map_candidates["sample"])
gene_col   = find_col(gc_df.columns, col_map_candidates["gene"])
gc_col     = find_col(gc_df.columns, col_map_candidates["gc"])

if sample_col is None or gene_col is None or gc_col is None:
    raise ValueError(f"Missing required columns in {gc_file}")

gc_df = gc_df.rename(columns={
    sample_col: "sample",
    gene_col: "gene_name",
    gc_col: "gc_percent"
})

gc_df["sample"] = gc_df["sample"].astype(str).str.strip()
gc_df["gene_name"] = gc_df["gene_name"].astype(str).str.strip()
gc_df["gc_percent"] = pd.to_numeric(gc_df["gc_percent"], errors="coerce")

# Convert if needed
if gc_df["gc_percent"].max(skipna=True) <= 1.05:
    gc_df["gc_percent"] *= 100

if REMOVE_RRNA_TRNA:
    gc_df = gc_df[
        ~gc_df["gene_name"].str.lower().str.contains(r"trna|rrna|rna", regex=True)
    ].copy()

# RAW GC (no normalization)
gc_df["plot_value"] = gc_df["gc_percent"]
gc_ylabel = "GC content (%)"

gc_df["group_key"] = gc_df["sample"].map(abbrev)
gc_df["species"] = gc_df["group_key"].map(sample_name_map).fillna(gc_df["sample"])

gc_subset = gc_df[gc_df["group_key"].isin(set(aurelia) | set(others))].copy()
gc_aurelia_df = gc_subset[gc_subset["group_key"].isin(aurelia)].copy()
gc_others_df  = gc_subset[gc_subset["group_key"].isin(others)].copy()

# =========================================================
# GENE ORDER
# =========================================================
all_genes = sorted(
    set(len_df["gene_name"]).intersection(set(gc_df["gene_name"]))
)

x = np.arange(len(all_genes))
rng = np.random.default_rng(3)

# =========================================================
# FIGURE
# =========================================================
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(max(18, 0.65 * len(all_genes)), 15),
    sharex=True
)

# =========================================================
# PANEL 1: GENE LENGTH
# =========================================================
len_box_data = [
    len_aurelia_df.loc[len_aurelia_df["gene_name"] == g, "plot_length"].dropna().values
    for g in all_genes
]

ax1.boxplot(
    len_box_data,
    positions=x,
    widths=0.6,
    showfliers=False,
    patch_artist=True,
    boxprops=dict(facecolor="lightgray", edgecolor="black"),
    medianprops=dict(color="black"),
)

# scatter
for i, g in enumerate(all_genes):
    vals = len_aurelia_df.loc[len_aurelia_df["gene_name"] == g, "plot_length"].dropna().values
    if len(vals):
        ax1.scatter(np.full(len(vals), i) + rng.uniform(-0.15,0.15,len(vals)),
                    vals, color="gray", s=22, alpha=0.7)

for sp, col in outgroup_colors.items():
    sub = len_others_df[len_others_df["species"] == sp]
    for i, g in enumerate(all_genes):
        vals = sub.loc[sub["gene_name"] == g, "plot_length"].dropna().values
        if len(vals):
            ax1.scatter(np.full(len(vals), i)+rng.uniform(-0.15,0.15,len(vals)),
                        vals, color=col, s=60, edgecolor="black")

ax1.axhline(1, linestyle="--", color="black")
ax1.set_ylabel(len_ylabel, fontsize=32)

# =========================================================
# PANEL 2: GC
# =========================================================
gc_box_data = [
    gc_aurelia_df.loc[gc_aurelia_df["gene_name"] == g, "plot_value"].dropna().values
    for g in all_genes
]

ax2.boxplot(gc_box_data, positions=x, widths=0.6, patch_artist=True,
            boxprops=dict(facecolor="lightgray", edgecolor="black"))

# scatter
for i, g in enumerate(all_genes):
    vals = gc_aurelia_df.loc[gc_aurelia_df["gene_name"] == g, "plot_value"].dropna().values
    if len(vals):
        ax2.scatter(np.full(len(vals), i)+rng.uniform(-0.15,0.15,len(vals)),
                    vals, color="gray", s=22)

for sp, col in outgroup_colors.items():
    sub = gc_others_df[gc_others_df["species"] == sp]
    for i, g in enumerate(all_genes):
        vals = sub.loc[sub["gene_name"] == g, "plot_value"].dropna().values
        if len(vals):
            ax2.scatter(np.full(len(vals), i)+rng.uniform(-0.15,0.15,len(vals)),
                        vals, color=col, s=55, edgecolor="black")

ax2.set_ylabel(gc_ylabel, fontsize=32)
ax2.set_xlabel("Gene name", fontsize=32)

# improved x labels
ax2.set_xticks(x)
ax2.set_xticklabels(all_genes, rotation=90, fontsize=24, color="black")

# =========================================================
# LEGEND (FINAL FIX)
# =========================================================
handles = [Line2D([0],[0], marker="o", color="gray", linestyle="None", markersize=9, label="Aurelia")]

for sp, col in outgroup_colors.items():
    handles.append(Line2D([0],[0], marker="o", color=col,
                          markeredgecolor="black", linestyle="None",
                          markersize=9, label=sp))

legend = ax1.legend(
    handles=handles,
    bbox_to_anchor=(1.02, 1.0),
    loc="upper left",
    fontsize=26,          # 🔥 bigger font
    title_fontsize=24,    # if you add title later
    borderpad=1.5,        # more padding inside box
    labelspacing=1.5,     # space between items
    handletextpad=1.0,    # space between marker and text
    markerscale=1.5       # bigger legend markers
)

for text in legend.get_texts():
    text.set_style("italic")

legend.get_frame().set_edgecolor("black")
legend.get_frame().set_linewidth(1.5)
legend.get_frame().set_facecolor("white")


# =========================================================
# FINAL
# =========================================================
plt.tight_layout()

# ✅ increase gap between the two plots
plt.subplots_adjust(
    hspace=0.25,   # 🔥 increase this (0.20–0.35 ideal)
    right=0.80
)

plt.savefig(OUT_FIG, dpi=600, bbox_inches="tight")
plt.close()

print("Saved:", OUT_FIG)

print("Saved:", OUT_FIG)