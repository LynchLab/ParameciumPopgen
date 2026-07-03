#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ============================================================
# Global style
# ============================================================
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["text.color"] = "#000000"
plt.rcParams["axes.labelcolor"] = "#000000"
plt.rcParams["xtick.color"] = "#000000"
plt.rcParams["ytick.color"] = "#000000"
plt.rcParams["axes.edgecolor"] = "#000000"

TITLE_FONTSIZE = 21
AXIS_LABEL_FONTSIZE = 20
TICK_FONTSIZE_X = 17
TICK_FONTSIZE_Y = 17
LEGEND_FONTSIZE = 17

# ============================================================
# Shared colors
# ============================================================
outgroup_colors = {
    "P. caudatum": "#e69138",
    "P. multimicronucleatum": "#9874f4",
    "P. fokini": "#ffe599",
    "P. bursaria": "#cc0000",
    "P. tribursaria": "#88c96b",
}
aurelia_point_color = "gray"
aurelia_box_color = "lightgray"

# ============================================================
# Figure layout
# ============================================================
fig = plt.figure(figsize=(18, 21))
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.55, 1.0],
    width_ratios=[1, 1],
    hspace=0.58,
    wspace=0.24
)

ax_top = fig.add_subplot(gs[0, :])
ax_bl = fig.add_subplot(gs[1, 0])
ax_br = fig.add_subplot(gs[1, 1])

# ============================================================
# PANEL A: Gene presence/absence heatmap
# ============================================================
INPUT_FILE = "gene_presence_absence_matrix_2.csv"

heat_df = pd.read_csv(INPUT_FILE, index_col=0)

if "orfX" in heat_df.index:
    heat_df = heat_df.drop(index="orfX")

drop_cols = [
    c for c in heat_df.columns
    if ("tetrahymena" in c.lower())
    or ("gigas" in c.lower())
    or ("pquin" in c.lower())
    or ("calkinski" in c.lower())
]
if drop_cols:
    heat_df = heat_df.drop(columns=drop_cols)

def to_bin(x):
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return 1
    elif s in {"0", "false", "f", "no", "n", ""}:
        return 0
    else:
        try:
            return 1 if float(s) != 0 else 0
        except Exception:
            return 0

# fixed deprecation warning
heat_df = heat_df.apply(lambda col: col.map(to_bin))

species_fullname_map = {
    "quindecaurelia": "P. quindecaurelia",
    "quadecaurelia": "P. quadecaurelia",
    "tredecaurelia": "P. tredecaurelia",
    "dodecaurelia": "P. dodecaurelia",
    "decaurelia": "P. decaurelia",
    "novaurelia": "P. novaurelia",
    "octaurelia": "P. octaurelia",
    "sexaurelia": "P. sexaurelia",
    "pentaurelia": "P. pentaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "biaurelia": "P. biaurelia",
    "primaurelia": "P. primaurelia",
    "sonneborni": "P. sonneborni",
    "jenningsi": "P. jenningsi",
    "Psche": "P. schewiakoffi",
    "caudatum": "P. caudatum",
    "multimicronucleatum": "P. multimicronucleatum",
    "multi": "P. multimicronucleatum",
    "fokini": "P. fokini",
    "bursaria": "P. bursaria",
    "tribursaria": "P. tribursaria",
}

compiled_species = [
    (k, re.compile(rf"(?i)(?:^|[^a-z]){re.escape(k)}(?:$|[^a-z])"))
    for k in sorted(species_fullname_map.keys(), key=len, reverse=True)
]

def rename_species(col):
    low = col.lower()
    for k, rx in compiled_species:
        if rx.search(low):
            return species_fullname_map[k]
    return col

heat_df.columns = [rename_species(c) for c in heat_df.columns]

aurelia_order = [
    "P. primaurelia", "P. biaurelia", "P. triaurelia", "P. tetraurelia",
    "P. pentaurelia", "P. sexaurelia", "P. octaurelia", "P. novaurelia",
    "P. decaurelia", "P. dodecaurelia", "P. tredecaurelia",
    "P. quadecaurelia", "P. quindecaurelia", "P. sonneborni",
    "P. jenningsi", "P. schewiakoffi"
]
group_order = [
    aurelia_order,
    ["P. caudatum"],
    ["P. multimicronucleatum"],
    ["P. fokini"],
    ["P. bursaria", "P. tribursaria"],
]

ordered_cols = []
seen = set()
for group in group_order:
    for sp in group:
        if sp in heat_df.columns and sp not in seen:
            ordered_cols.append(sp)
            seen.add(sp)

ordered_cols += [c for c in heat_df.columns if c not in seen]
heat_df = heat_df[ordered_cols]

gene_order = [
    "atp9", "cob", "cox1", "cox2", "nadh1-a", "nadh1-b", "nadh2", "nadh3",
    "nadh4", "nadh4L", "nadh5", "nadh6", "nadh7", "nadh9", "nadh10",
    "rpL2", "rpL6", "rpL14", "rpL16", "rpL32", "rpS3", "rpS12", "rpS13",
    "rpS14", "rpS19", "yejR", "ymf56", "ymf57", "ymf59", "ymf61", "ymf63",
    "ymf64", "ymf65", "ymf66", "ymf67", "ymf68", "ymf76", "ymf78", "ymf79",
    "ymf80", "ymf83", "ymf84", "orfL", "orfY", "orfZ"
]
final_rows = [g for g in gene_order if g in heat_df.index] + [g for g in heat_df.index if g not in gene_order]
heat_df = heat_df.loc[final_rows]

heat_cmap = sns.color_palette(["#f4e04d", "#4c78a8"])
sns.heatmap(
    heat_df,
    cmap=heat_cmap,
    cbar=False,
    linewidths=0.8,
    linecolor="black",
    xticklabels=True,
    yticklabels=True,
    ax=ax_top
)

ax_top.set_xlabel("Species", fontsize=AXIS_LABEL_FONTSIZE + 1, color="#000000", labelpad=16)
ax_top.set_ylabel("Genes", fontsize=AXIS_LABEL_FONTSIZE + 1, color="#000000", labelpad=12)

ax_top.tick_params(
    axis="x",
    rotation=45,
    labelsize=TICK_FONTSIZE_X,
    colors="#000000",
    pad=8
)
ax_top.tick_params(
    axis="y",
    labelsize=TICK_FONTSIZE_Y + 1,
    colors="#000000",
    pad=4
)

for tick in ax_top.get_xticklabels():
    tick.set_horizontalalignment("right")

# heat_legend = [
#     Patch(facecolor="#f4e04d", edgecolor="black", label="Absent"),
#     Patch(facecolor="#4c78a8", edgecolor="black", label="Present"),
# ]
# leg_a = ax_top.legend(
#     handles=heat_legend,
#     loc="upper right",
#     bbox_to_anchor=(1.18, 1.04),
#     frameon=True,
#     fontsize=LEGEND_FONTSIZE + 3,
#     borderpad=1.2,
#     handlelength=2.0,
#     handleheight=1.5,
#     labelspacing=0.8
# )

# frame = leg_a.get_frame()
# frame.set_edgecolor("black")
# frame.set_linewidth(1.5)
# frame.set_alpha(1.0)

# ============================================================
# Shared input for PANEL B
# ============================================================
coords_file = "combined_gene_coords_table.csv"
coords_df = pd.read_csv(coords_file)
coords_df["sample"] = coords_df["sample"].astype(str).str.strip()
coords_df["gene_name"] = coords_df["gene_name"].astype(str).str.strip()
coords_df = coords_df[~coords_df["gene_name"].str.lower().str.contains(r"trna|rrna|rna", regex=True)].copy()
coords_df["normalized_length"] = coords_df.groupby("gene_name")["length"].transform(lambda x: x / x.median())

aurelia_codes = [
    "tetr", "tria", "pent", "sexa", "octa", "dode", "tred",
    "biau", "deca", "jenn", "nova", "prima", "psch", "quad", "sonn"
]
other_codes = ["caud", "mult", "foki", "burs", "tbur", "calk"]

sample_name_map = {
    "caud": "P. caudatum",
    "burs": "P. bursaria",
    "tbur": "P. tribursaria",
    "mult": "P. multimicronucleatum",
    "foki": "P. fokini",
    "calk": "P. calkinski",
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
    "psch": "P. schewiakoffi",
    "quad": "P. quadecaurelia",
    "sonn": "P. sonneborni",
}
coords_df["species"] = coords_df["sample"].map(sample_name_map).fillna(coords_df["sample"])

aurelia_df = coords_df[coords_df["sample"].isin(aurelia_codes)].copy()
others_df = coords_df[coords_df["sample"].isin(other_codes)].copy()

# ============================================================
# PANEL B: Gene length differences
# ============================================================
selected_genes = ["cox2", "cob", "nadh1-b", "nadh10", "nadh5", "rpS19", "ymf66"]
x = np.arange(len(selected_genes))

box_data = [
    aurelia_df.loc[aurelia_df["gene_name"] == g, "normalized_length"].dropna().values
    for g in selected_genes
]

ax_bl.boxplot(
    box_data,
    positions=x,
    widths=0.6,
    showfliers=False,
    patch_artist=True,
    boxprops=dict(facecolor=aurelia_box_color, edgecolor="black", linewidth=0.8),
    medianprops=dict(color="black", linewidth=1.1),
    whiskerprops=dict(color="black", linewidth=0.8),
    capprops=dict(color="black", linewidth=0.8),
)

rng = np.random.default_rng(3)

for i, g in enumerate(selected_genes):
    vals = aurelia_df.loc[aurelia_df["gene_name"] == g, "normalized_length"].dropna().values
    if len(vals):
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax_bl.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            color=aurelia_point_color,
            s=18,
            alpha=0.7,
            zorder=3
        )

for sp, col in outgroup_colors.items():
    sub_sp = others_df[others_df["species"] == sp]
    for i, g in enumerate(selected_genes):
        vals = sub_sp.loc[sub_sp["gene_name"] == g, "normalized_length"].dropna().values
        if len(vals):
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax_bl.scatter(
                np.full(len(vals), i) + jitter,
                vals,
                color=col,
                s=49,
                edgecolor="black",
                linewidth=0.4,
                zorder=4
            )

ax_bl.axhline(1.0, linestyle="--", linewidth=1, color="black", alpha=0.6)
ax_bl.set_xticks(x)
ax_bl.set_xticklabels(
    selected_genes,
    rotation=90,
    fontsize=TICK_FONTSIZE_X + 5,

    color="#000000"
)
ax_bl.tick_params(axis="x", pad=10, length=6, width=1.2)
ax_bl.set_ylabel("Normalized gene length", fontsize=AXIS_LABEL_FONTSIZE + 3, color="#000000")
ax_bl.set_xlabel("")
ax_bl.set_ylim(0.9, 1.7)
ax_bl.tick_params(axis="y", labelsize=TICK_FONTSIZE_Y + 1, colors="#000000")
ax_bl.margins(x=0.05)

# ============================================================
# PANEL C: Median GC compartments
# ============================================================
gc_data = [
    {"Species": "P. biaurelia",           "Protein-coding": 32.83, "rRNA": 35.86, "tRNA": 41.36, "Intergenic": 32.49},
    {"Species": "P. decaurelia",          "Protein-coding": 38.91, "rRNA": 36.37, "tRNA": 43.73, "Intergenic": 40.84},
    {"Species": "P. dodecaurelia",        "Protein-coding": 40.30, "rRNA": 36.52, "tRNA": 43.73, "Intergenic": 40.58},
    {"Species": "P. jenningsi",           "Protein-coding": 39.76, "rRNA": 36.80, "tRNA": 44.67, "Intergenic": 39.87},
    {"Species": "P. novaurelia",          "Protein-coding": 40.08, "rRNA": 35.97, "tRNA": 43.73, "Intergenic": 40.14},
    {"Species": "P. octaurelia",          "Protein-coding": 39.69, "rRNA": 35.77, "tRNA": 43.73, "Intergenic": 39.24},
    {"Species": "P. pentaurelia",         "Protein-coding": 39.63, "rRNA": 35.88, "tRNA": 43.73, "Intergenic": 39.85},
    {"Species": "P. quindecaurelia",      "Protein-coding": 33.34, "rRNA": 35.60, "tRNA": 44.44, "Intergenic": 29.07},
    {"Species": "P. primaurelia",         "Protein-coding": 35.85, "rRNA": 35.63, "tRNA": 40.84, "Intergenic": 40.57},
    {"Species": "P. schewiakoffi",        "Protein-coding": 37.61, "rRNA": 36.39, "tRNA": 43.73, "Intergenic": 37.49},
    {"Species": "P. quadecaurelia",       "Protein-coding": 38.71, "rRNA": 36.10, "tRNA": 38.18, "Intergenic": 31.77},
    {"Species": "P. sexaurelia",          "Protein-coding": 37.68, "rRNA": 35.30, "tRNA": 39.09, "Intergenic": 36.94},
    {"Species": "P. sonneborni",          "Protein-coding": 40.02, "rRNA": 36.43, "tRNA": 43.73, "Intergenic": 40.76},
    {"Species": "P. tetraurelia",         "Protein-coding": 41.54, "rRNA": 35.98, "tRNA": 41.36, "Intergenic": 38.44},
    {"Species": "P. triaurelia",          "Protein-coding": 35.65, "rRNA": 35.58, "tRNA": 40.84, "Intergenic": 33.66},
    {"Species": "P. caudatum",            "Protein-coding": 20.53, "rRNA": 32.03, "tRNA": 40.84, "Intergenic": 15.49},
    {"Species": "P. multimicronucleatum", "Protein-coding": 19.15, "rRNA": 31.26, "tRNA": 40.84, "Intergenic": 13.87},
    {"Species": "P. fokini",              "Protein-coding": 20.80, "rRNA": 31.07, "tRNA": 39.50, "Intergenic": 12.52},
    {"Species": "P. bursaria",            "Protein-coding": 14.70, "rRNA": 26.24, "tRNA": 37.70, "Intergenic": 10.35},
    {"Species": "P. tribursaria",         "Protein-coding": 14.55, "rRNA": 26.94, "tRNA": 38.19, "Intergenic": 8.57},
]

gc_df = pd.DataFrame(gc_data)
outgroups = list(outgroup_colors.keys())
gc_df["Group"] = gc_df["Species"].apply(lambda s: "Outgroup" if s in outgroups else "Aurelia")

gc_long = gc_df.melt(
    id_vars=["Species", "Group"],
    value_vars=["Protein-coding", "Intergenic", "rRNA", "tRNA"],
    var_name="Region",
    value_name="GC"
).dropna()

region_order = ["Protein-coding", "Intergenic", "rRNA", "tRNA"]

# custom spacing with extra gap between first two categories
custom_positions = {
    "Protein-coding": 0.0,
    "Intergenic": 1.8,
    "rRNA": 3.0,
    "tRNA": 4.2
}
gc_long["xpos"] = gc_long["Region"].map(custom_positions)

# Aurelia boxplots by manual positions
box_data_gc = [
    gc_long.loc[(gc_long["Group"] == "Aurelia") & (gc_long["Region"] == r), "GC"].dropna().values
    for r in region_order
]
positions_gc = [custom_positions[r] for r in region_order]

ax_br.boxplot(
    box_data_gc,
    positions=positions_gc,
    widths=0.55,
    showfliers=False,
    patch_artist=True,
    boxprops=dict(facecolor=aurelia_box_color, edgecolor="black", linewidth=0.8),
    medianprops=dict(color="black", linewidth=1.1),
    whiskerprops=dict(color="black", linewidth=0.8),
    capprops=dict(color="black", linewidth=0.8),
    zorder=1
)

# Aurelia points
for r in region_order:
    vals = gc_long.loc[(gc_long["Group"] == "Aurelia") & (gc_long["Region"] == r), "GC"].dropna().values
    if len(vals):
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax_br.scatter(
            np.full(len(vals), custom_positions[r]) + jitter,
            vals,
            color=aurelia_point_color,
            s=18,
            alpha=0.7,
            zorder=2
        )

# Outgroup points
for sp, col in outgroup_colors.items():
    sub_sp = gc_long[gc_long["Species"] == sp]
    for r in region_order:
        vals = sub_sp.loc[sub_sp["Region"] == r, "GC"].dropna().values
        if len(vals):
            jitter = rng.uniform(-0.12, 0.12, size=len(vals))
            ax_br.scatter(
                np.full(len(vals), custom_positions[r]) + jitter,
                vals,
                color=col,
                s=49,
                edgecolor="black",
                linewidth=0.4,
                zorder=3
            )

ax_br.set_xlabel("")
ax_br.set_ylabel("Median GC content (%)", fontsize=AXIS_LABEL_FONTSIZE + 3, color="#000000")
ax_br.set_ylim(0, 50)
ax_br.set_xlim(-0.6, 4.8)
ax_br.set_xticks(positions_gc)
ax_br.set_xticklabels(
    region_order,
    fontsize=TICK_FONTSIZE_X + 4,

    color="#000000"
)

for label in ax_br.get_xticklabels():
    label.set_rotation(20)
    label.set_horizontalalignment("right")

ax_br.tick_params(
    axis="x",
    colors="#000000",
    pad=12,
    length=6,
    width=1.2
)
ax_br.tick_params(axis="y", labelsize=TICK_FONTSIZE_Y + 2, colors="#000000")
ax_br.margins(x=0.10)

# ============================================================
# Shared legend for B and C
# ============================================================
shared_handles = [
    Line2D([0], [0], marker="o", color=aurelia_point_color, linestyle="None", markersize=6, label="Aurelia"),
]
for sp, col in outgroup_colors.items():
    shared_handles.append(
        Line2D([0], [0], marker="o", color=col, markeredgecolor="black",
               linestyle="None", markersize=6, label=sp)
    )

# leg_bc = fig.legend(
#     handles=shared_handles,
#     loc="lower center",
#     bbox_to_anchor=(0.74, 0.02),
#     ncol=3,
#     frameon=True,
#     fontsize=LEGEND_FONTSIZE + 1,
#     borderpad=0.4,
#     handlelength=1.0,
#     handleheight=0.8,
#     labelspacing=0.3,
#     columnspacing=0.8
# )
#
# frame = leg_bc.get_frame()
# frame.set_edgecolor("black")
# frame.set_linewidth(1.4)
# frame.set_alpha(1.0)

# for txt in leg_bc.get_texts():
#     txt.set_fontstyle("italic")
#     txt.set_color("black")

# ============================================================
# Final spacing and save
# ============================================================
plt.subplots_adjust(top=0.96, bottom=0.11, left=0.07, right=0.98)
plt.savefig("combined_three_panel_figure.png", dpi=600, bbox_inches="tight")
plt.savefig("combined_three_panel_figure.pdf", bbox_inches="tight")
plt.show()