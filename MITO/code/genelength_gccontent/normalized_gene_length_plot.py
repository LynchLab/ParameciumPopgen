import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------- Load data ----------
file_path = "combined_gene_coords_table.csv"
df = pd.read_csv(file_path)

# Clean names
df["sample"] = df["sample"].astype(str).str.strip()
df["gene_name"] = df["gene_name"].astype(str).str.strip()

# ---------- Remove tRNA / rRNA ----------
mask = ~df["gene_name"].str.lower().str.contains("trna|rrna|rna")
df = df[mask].copy()

# ---------- Normalize each gene by its median ----------
df["normalized_length"] = df.groupby("gene_name")["length"].transform(
    lambda x: x / x.median()
)

# ---------- Define species groups ----------
aurelia = [
    "tetr", "tria", "pent", "sexa", "octa", "dode", "tred",
    "biau", "deca", "jenn", "nova", "prima", "psch", "quad", "sonn"
]
others = ["caud", "mult", "foki", "burs", "tbur", "calk"]

# ---------- Map sample codes to proper species names ----------
sample_name_map = {
    "caud": "P. caudatum",
    "burs": "P. bursaria",
    "tbur": "P. tribursaria",
    "mult": "P. multimicronucleatum",
    "foki": "P. fokini",
    "calk": "P. calkinsi",

    # Aurelia (not shown individually in legend)
    "tetr": "P. tetraurelia",
    "tria": "P. triaurelia",
    "pent": "P. pentaurelia",
    "sexa": "P. sexaurelia",
    "octa": "P. octaurelia",
    "dode": "P. dodecaurelia",
    "tred": "P. tredecaurelia",
    "biau": "P. biaurelia",
    "deca": "P. decaurelia",
    "jenn": "P. jenningsi",
    "nova": "P. novaurelia",
    "prima": "P. primaurelia",
    "psch": "P. pschorrigi",
    "quad": "P. quadecaurelia",
    "sonn": "P. sonneborni"
}

df["species"] = df["sample"].map(sample_name_map).fillna(df["sample"])

aurelia_df = df[df["sample"].isin(aurelia)].copy()
others_df = df[df["sample"].isin(others)].copy()

# ---------- Gene order ----------
gene_order = sorted(df["gene_name"].unique())

# ---------- Plot ----------
plt.figure(figsize=(22, 8))

# Aurelia boxplots (gray)
sns.boxplot(
    data=aurelia_df,
    x="gene_name",
    y="normalized_length",
    order=gene_order,
    color="lightgray",
    showfliers=False,
    zorder=1
)

# Aurelia points (gray)
sns.stripplot(
    data=aurelia_df,
    x="gene_name",
    y="normalized_length",
    order=gene_order,
    color="gray",
    size=4,
    jitter=0.18,
    alpha=0.7,
    zorder=2
)

# Non-Aurelia colored points
palette = sns.color_palette("tab10", n_colors=len(others))
sns.stripplot(
    data=others_df,
    x="gene_name",
    y="normalized_length",
    order=gene_order,
    hue="species",
    hue_order=[sample_name_map[s] for s in others],
    dodge=True,
    size=5,
    jitter=0.18,
    palette=palette,
    zorder=3
)

# ---------- Styling ----------
plt.title(
    "Relative Mitochondrial Gene Lengths\nAurelia (gray) vs Non-Aurelia (colored)",
    fontsize=18,
    weight="bold"
)
plt.xlabel("Gene", fontsize=13)
plt.ylabel("Relative gene length (length / median per gene)", fontsize=13)
plt.xticks(rotation=90, fontsize=11)
plt.yticks(fontsize=11)

# Reference line at median
plt.axhline(1.0, linestyle="--", color="black", linewidth=1, alpha=0.6)

# ---------- Custom legend ----------
handles, labels = plt.gca().get_legend_handles_labels()

aurelia_handle = Line2D(
    [0], [0],
    marker="o",
    color="gray",
    linestyle="None",
    markersize=7,
    label="Aurelia"
)

handles = [aurelia_handle] + handles
labels = ["Aurelia"] + labels

plt.legend(
    handles,
    labels,
    title="Species",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0.
)

plt.tight_layout()

# ---------- Save (PNG, high resolution) ----------
out_png = "aurelia_normalized_gray_boxplots_colored_others.png"
plt.savefig(out_png, dpi=600)
plt.close()

print(f"✅ Saved: {out_png} (PNG, 600 dpi, publication-ready)")
