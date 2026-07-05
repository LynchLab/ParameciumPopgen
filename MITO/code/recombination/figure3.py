import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import os


# -----------------------------
# Input
# -----------------------------
phi_path = "phi_full_summary.csv"
print(os.path.exists(phi_path))
# Output
out_png = "FigS_recombination_tests_all_side_by_side.png"
out_pdf = "FigS_recombination_tests_all_side_by_side.pdf"

# -----------------------------
# Load
# -----------------------------
phi = pd.read_csv(phi_path)

# -----------------------------
# Species labels
# -----------------------------
species_labels = {
    "biaurelia": "P. biaurelia",
    "bursaria": "P. bursaria",
    "caudatum": "P. caudatum",
    "multimicronucleatum": "P. multimicronucleatum",
    "multi": "P. multimicronucleatum",
    "fokini": "P. fokini",
    "quindecaurelia": "P. quindecaurelia",
    "pquin": "P. quindecaurelia",
    "primaurelia": "P. primaurelia",
    "sexaurelia": "P. sexaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "tribursaria": "P. tribursaria",
    # ✅ ADD THESE
    "tredecaurelia": "P. tredecaurelia",
    "quadecaurelia": "P. quadecaurelia",
}

# -----------------------------
# Helpers
# -----------------------------
def norm_species(s):
    s = str(s).strip()
    s = s.replace("Paramecium", "").strip()
    s = s.replace("P.", "").strip()
    s = s.replace(" ", "").replace(".", "")
    return s.lower()

def gene_group_rank(gene):
    g = str(gene).lower()
    if g.startswith("cox"):
        return (0, g)
    if g == "cob" or g.startswith("cytb"):
        return (1, g)
    if g.startswith("nadh"):
        return (2, g)
    if g.startswith("rpl"):
        return (3, g)
    if g.startswith("rps"):
        return (4, g)
    if g.startswith("orf"):
        return (5, g)
    if g.startswith("ymf"):
        return (6, g)
    return (7, g)

def find_col(df, must_contain, optional_contain=None):
    must = [x.lower() for x in (must_contain if isinstance(must_contain, (list, tuple)) else [must_contain])]
    opt = [x.lower() for x in (optional_contain if optional_contain else [])]
    for c in df.columns:
        cl = c.lower()
        if all(m in cl for m in must) and all(o in cl for o in opt):
            return c
    return None

def to_numeric_clean(series):
    s = series.astype(str).str.strip()
    s = s.replace({
        "NA": np.nan, "na": np.nan, "NaN": np.nan, "nan": np.nan,
        "None": np.nan, "none": np.nan, "": np.nan
    })
    s = s.str.replace(",", "", regex=False)
    return pd.to_numeric(s, errors="coerce")

# -----------------------------
# Global style
# -----------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
})

# -----------------------------
# Validate required columns
# -----------------------------
required = {"Species", "Gene"}
missing_req = [c for c in required if c not in phi.columns]
if missing_req:
    raise KeyError(
        f"Missing required columns in phi_full_summary.csv: {missing_req}\n"
        f"Columns: {list(phi.columns)}"
    )

# -----------------------------
# Standardize keys
# -----------------------------
phi2 = phi.copy()
phi2["species_key"] = phi2["Species"].apply(norm_species)
phi2["gene"] = phi2["Gene"].astype(str)

# -----------------------------
# Species order
# -----------------------------
preferred_species = [
    "biaurelia", "bursaria", "caudatum", "fokini", "primaurelia",
    "sexaurelia", "tetraurelia", "triaurelia", "tribursaria", "multimicronucleatum"
]
present_species = list(pd.unique(phi2["species_key"]))
species_order = [s for s in preferred_species if s in present_species] + \
                [s for s in present_species if s not in preferred_species]

# -----------------------------
# Gene order
# -----------------------------
gene_order = sorted(pd.unique(phi2["gene"]), key=gene_group_rank)

# -----------------------------
# Detect p-value columns
# -----------------------------
col_nss = find_col(phi2, "nss")
col_maxchi2 = find_col(phi2, ["maxchi"], ["2"]) or find_col(phi2, "maxchi2") or find_col(phi2, "maxchi")
col_phiperm = find_col(phi2, ["phi"], ["perm"]) or find_col(phi2, "phi_perm")
col_phinorm = find_col(phi2, ["phi"], ["norm"]) or find_col(phi2, "phi_norm")

test_cols = [
    ("NSS", col_nss),
    ("MaxChi2", col_maxchi2),
    ("PHI_perm", col_phiperm),
    ("PHI_norm", col_phinorm),
]

missing = [name for name, c in test_cols if c is None]
if missing:
    raise KeyError(
        f"Could not find columns for: {missing}\n"
        f"Available columns: {list(phi2.columns)}\n\n"
        "If your headers differ, rename them or edit the find_col() logic."
    )

# -----------------------------
# Convert p-value columns to numeric
# -----------------------------
for _, pcol in test_cols:
    phi2[pcol] = to_numeric_clean(phi2[pcol])

# -----------------------------
# Build categorical matrices
# 0 = no data, 1 = absent, 2 = present
# -----------------------------
sig_thr = 0.05
mats = {}

for test_name, pcol in test_cols:
    piv = phi2.pivot_table(index="gene", columns="species_key", values=pcol, aggfunc="min")
    piv = piv.reindex(index=gene_order, columns=species_order)

    cat = pd.DataFrame(0, index=piv.index, columns=piv.columns, dtype=int)
    cat[(piv >= sig_thr) & (~piv.isna())] = 1
    cat[(piv < sig_thr) & (~piv.isna())] = 2

    mats[test_name] = cat

# -----------------------------
# Plot
# -----------------------------
cmap = ListedColormap(["#d9d9d9", "#f4a742", "#2ca25f"])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

fig, axes = plt.subplots(
    1, 4,
    figsize=(24, 11),
    gridspec_kw={"wspace": 0.14}
)

for ax, (test_name, _) in zip(axes, test_cols):
    cat = mats[test_name]

    ax.imshow(
        cat.values,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest"
    )

    # X labels: species
    ax.set_xticks(range(len(cat.columns)))
    ax.set_xticklabels(
        [species_labels.get(s, s) for s in cat.columns],
        rotation=45,
        ha="right",
        fontsize=12,
        color="black"
    )

    # Y labels: genes only on first panel
    ax.set_yticks(range(len(cat.index)))
    if ax is axes[0]:
        ax.set_yticklabels(
            cat.index,
            fontsize=14,
            color="black"
        )
        ax.set_ylabel("Gene", fontsize=20, color="black", labelpad=10)
    else:
        ax.set_yticklabels([])

    # Panel title
    ax.set_title(
        test_name,
        fontsize=20,
        color="black",
        pad=12
    )

    # Tick appearance
    ax.tick_params(axis="x", colors="black", labelsize=12)
    ax.tick_params(axis="y", colors="black", labelsize=12)

    # Black borders between cells
    ax.set_xticks(np.arange(-0.5, cat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, cat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=0.4)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Axes border
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)

# -----------------------------
# Legend on the side with box
# -----------------------------
legend_handles = [
    Patch(facecolor="#2ca25f", edgecolor="black",
          label="Nominal recombination signal (P < 0.05)"),
    Patch(facecolor="#f4a742", edgecolor="black",
          label="No significant signal (P ≥ 0.05)"),
    Patch(facecolor="#d9d9d9", edgecolor="black",
          label="Insufficient data"),
]

legend = fig.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=True,
    fontsize=13,
    borderpad=1.2,
    labelspacing=1.1,
    handletextpad=0.8
)

for text in legend.get_texts():
    text.set_color("black")
    text.set_fontstyle("normal")
    text.set_fontweight("normal")

legend.get_frame().set_edgecolor("black")
legend.get_frame().set_linewidth(1.2)
legend.get_frame().set_alpha(1.0)

# -----------------------------
# Final layout
# -----------------------------
plt.tight_layout()
plt.subplots_adjust(right=0.84)

# Save
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("Wrote:", out_png)
print("Wrote:", out_pdf)
print("Detected columns:", {name: col for name, col in test_cols})

# FigX_recombination_main.py
# MAIN FIGURE:
# (A) PHI permutation test (Pairwise Homoplasy Index; PHI) heatmap (gene × species)
#     - Heatmap colors: green=present, orange=absent, light grey=no data
#     - Black borders between heatmap boxes
#     - Genes ordered: cox* → cob → nadh* → rpS* → rpL* → orf* → ymf* (then anything else)
# (B) LDhat 4Ner (Hudson vs Wakeley) at smallest θ
#     - NO connecting lines between the two 4Ner estimates
#     - y-axis starts slightly below 0 so 0 points are visible
# - BOTH legends placed OUTSIDE the panels (right side)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# -----------------------------
# Input files
# -----------------------------
ldhat_path = "poster images/paper figures/data/recombination/pairwise_recombination_table_ldhat.csv"
phi_path   = "poster images/paper figures/data/recombination/phi_full_summary.csv"

# Output files
out_png = "poster images/paper figures/Images/FigX_PHIperm_LDhat_2methods_main.png"
out_pdf = "poster images/paper figures/Images/FigX_PHIperm_LDhat_2methods_main.pdf"

# -----------------------------
# Load data
# -----------------------------
ld = pd.read_csv(ldhat_path)
phi = pd.read_csv(phi_path)

# -----------------------------
# Species labels + colors
# -----------------------------
species_labels = {
    "biaurelia": "P. biaurelia",
    "bursaria": "P. bursaria",
    "caudatum": "P. caudatum",
    "multimicronucleatum": "P. multimicronucleatum",
    "multi": "P. multimicronucleatum",
    "fokini": "P. fokini",
    "quindecaurelia": "P. quindecaurelia",
    "pquin": "P. quindecaurelia",
    "primaurelia": "P. primaurelia",
    "sexaurelia": "P. sexaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "tribursaria": "P. tribursaria",
    "tredecaurelia": "P. tredecaurelia",
    "quadecaurelia": "P. quadecaurelia",
}

color_map = {
    "tribursaria": "#88c96b",
    "bursaria": "#cc0000",
    "multi": "#9874f4",
    "multimicronucleatum": "#9874f4",
    "fokini": "#ffe599",
    "caudatum": "#e69138",
    "sexaurelia": "#559738",
    "tetraurelia": "#0b5394",
    "quindecaurelia": "#d5ed8a",
    "pquin": "#d5ed8a",
    "biaurelia": "#ee79ac",
    "triaurelia": "#d5a6bd",
    "primaurelia": "#6fa8dc",
    "tredecaurelia": "#76a5af",
    "quadecaurelia": "#b4a7d6",
}

# -----------------------------
# Global style
# -----------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
})

# -----------------------------
# Helpers
# -----------------------------
def norm_species(s):
    s = str(s).strip()
    s = s.replace("Paramecium", "").strip()
    s = s.replace("P.", "").strip()
    s = s.replace(" ", "").replace(".", "")
    return s.lower()

def find_col(df, must_contain, optional_contain=None):
    must = [x.lower() for x in (must_contain if isinstance(must_contain, (list, tuple)) else [must_contain])]
    opt  = [x.lower() for x in (optional_contain if optional_contain else [])]
    for c in df.columns:
        cl = c.lower()
        if all(m in cl for m in must) and all(o in cl for o in opt):
            return c
    return None

def gene_group_rank(gene):
    g = str(gene).lower()
    if g.startswith("cox"):
        return (0, g)
    if g == "cob" or g.startswith("cytb"):
        return (1, g)
    if g.startswith("nadh"):
        return (2, g)
    if g.startswith("rps"):
        return (3, g)
    if g.startswith("rpl"):
        return (4, g)
    if g.startswith("orf"):
        return (5, g)
    if g.startswith("ymf"):
        return (6, g)
    return (7, g)

# -----------------------------
# Build species order
# -----------------------------
phi_species = set(phi["Species"].apply(norm_species)) if "Species" in phi.columns else set()
ld_species  = set(ld["species"].apply(norm_species)) if "species" in ld.columns else set()

preferred_species = [
    "biaurelia", "bursaria", "caudatum", "multimicronucleatum", "fokini",
    "quindecaurelia", "primaurelia", "sexaurelia", "tetraurelia",
    "triaurelia", "tribursaria", "tredecaurelia", "quadecaurelia"
]

species_order = [s for s in preferred_species if (s in phi_species) or (s in ld_species)]

# ============================================================
# PANEL A: PHI permutation heatmap
# ============================================================
phi2 = phi.copy()
phi2["species_key"] = phi2["Species"].apply(norm_species)
phi2["gene"] = phi2["Gene"].astype(str)

p_col = "PHI_perm"
piv = phi2.pivot_table(index="gene", columns="species_key", values=p_col, aggfunc="min")
piv = piv.reindex(columns=species_order)

sig_thr = 0.05
cat = pd.DataFrame(0, index=piv.index, columns=piv.columns, dtype=int)
cat[(piv >= sig_thr) & (~piv.isna())] = 1
cat[(piv < sig_thr) & (~piv.isna())] = 2

# keep heatmap readable
has_any_signal = (cat == 2).sum(axis=1) > 0
has_enough_data = (cat != 0).sum(axis=1) >= max(3, int(0.25 * len(species_order)))
cat = cat.loc[has_any_signal | has_enough_data]

# order genes biologically
cat = cat.loc[sorted(cat.index, key=gene_group_rank)]

heat_cmap = ListedColormap(["#d9d9d9", "#f4a742", "#2ca25f"])
heat_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], heat_cmap.N)

# ============================================================
# PANEL B: LDhat 4Ner (Hudson vs Wakeley) at smallest theta
# ============================================================
ld2 = ld.copy()
ld2["species_key"] = ld2["species"].apply(norm_species)

# keep gene conversion model
ld2 = ld2[ld2["Model"].astype(str).str.contains("gene", case=False, na=False)].copy()

theta_col = find_col(ld2, "θ") or find_col(ld2, "theta") or "θ"
hudson_col = find_col(ld2, ["4ner"], ["hudson"]) or find_col(ld2, ["4ner"], ["finite"])
wakeley_col = find_col(ld2, ["4ner"], ["wakely"]) or find_col(ld2, ["4ner"], ["1997"])

if hudson_col is None or wakeley_col is None:
    raise KeyError(
        "Could not identify both 4Ner columns.\n"
        f"Columns seen: {list(ld2.columns)}"
    )

ld2 = ld2[
    ld2[theta_col].notna() &
    ld2[hudson_col].notna() &
    ld2[wakeley_col].notna()
].copy()

ld2 = ld2[ld2["species_key"].isin(species_order)].copy()

ld_min = (
    ld2.sort_values(["species_key", theta_col])
       .groupby("species_key", as_index=False)
       .first()
)
ld_min = ld_min.set_index("species_key").reindex(species_order).reset_index()

ld_min["label"] = ld_min["species_key"].map(species_labels).fillna(ld_min["species_key"])
ld_min["color"] = ld_min["species_key"].map(color_map).fillna("#999999")

y_h = ld_min[hudson_col].astype(float).values
y_w = ld_min[wakeley_col].astype(float).values

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(21, 8))
gs = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.5], wspace=0.30)
# -----------------------------
# Panel A: heatmap
# -----------------------------
axA = fig.add_subplot(gs[0, 0])
axA.imshow(cat.values, aspect="auto", cmap=heat_cmap, norm=heat_norm, interpolation="nearest")

axA.set_xticks(range(len(cat.columns)))
axA.set_xticklabels(
    [species_labels.get(sp, sp) for sp in cat.columns],
    rotation=45,
    ha="right",
    fontsize=20,
    color="black"
)

axA.set_yticks(range(len(cat.index)))
axA.set_yticklabels(
    cat.index,
    fontsize=5,
    color="black"
)

axA.set_xlabel("Species", fontsize=20, color="black", labelpad=14)
axA.set_ylabel("Gene", fontsize=20, color="black", labelpad=14)

axA.tick_params(axis="x", labelsize=20, colors="black")
axA.tick_params(axis="y", labelsize=10, colors="black")

# black borders between heatmap cells
axA.set_xticks(np.arange(-0.5, cat.shape[1], 1), minor=True)
axA.set_yticks(np.arange(-0.5, cat.shape[0], 1), minor=True)
axA.grid(which="minor", color="black", linestyle="-", linewidth=0.35)
axA.tick_params(which="minor", bottom=False, left=False)

# -----------------------------
# Panel B: LDhat points
# -----------------------------
axB = fig.add_subplot(gs[0, 1])

x = np.arange(len(ld_min))
dx = 0.14
x_h = x - dx
x_w = x + dx

axB.scatter(
    x_h, y_h,
    s=110,
    c=ld_min["color"],
    marker="o",
    edgecolor="black",
    linewidth=0.7,
    zorder=2
)
axB.scatter(
    x_w, y_w,
    s=110,
    c=ld_min["color"],
    marker="s",
    edgecolor="black",
    linewidth=0.7,
    zorder=2
)

axB.axhline(0, color="black", linestyle="--", linewidth=0.8)
axB.set_xlabel("Species", fontsize=20, color="black", labelpad=14)
axB.set_xticks(x)
axB.set_xticklabels(
    ld_min["label"],
    rotation=45,
    ha="right",
    fontsize=20,
    color="black"
)

axB.set_ylabel(
    r"Population-scaled recombination rate (4$N_e$r)",
    fontsize=20,
    color="black",
    labelpad=14
)

axB.tick_params(axis="x", labelsize=20, colors="black")
axB.tick_params(axis="y", labelsize=16, colors="black")

ymax = np.nanmax(np.concatenate([y_h, y_w])) if len(y_h) else 1.0
if not np.isfinite(ymax):
    ymax = 1.0
axB.set_ylim(-0.08, ymax * 1.10)

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# -----------------------------
# Legends outside on the right
# -----------------------------
heat_legend = [
    Patch(facecolor="#2ca25f", edgecolor="black",
          label="Nominal PHI signal (P < 0.05)"),
    Patch(facecolor="#f4a742", edgecolor="black",
          label="No significant PHI signal (P ≥ 0.05)"),
    Patch(facecolor="#d9d9d9", edgecolor="black",
          label="Insufficient data"),
]

method_legend = [
    Line2D([0], [0], marker="o", color="black", markerfacecolor="white",
           markersize=10, linestyle="None", label="Hudson estimator"),
    Line2D([0], [0], marker="s", color="black", markerfacecolor="white",
           markersize=10, linestyle="None", label="Wakeley estimator"),
]

plt.subplots_adjust(right=0.78)

legend1 = fig.legend(
    handles=heat_legend,
    frameon=True,
    loc="center left",
    bbox_to_anchor=(0.80, 0.72),
    fontsize=16,
    borderpad=1.6,
    labelspacing=1.4,
    handletextpad=0.8
)

legend2 = fig.legend(
    handles=method_legend,
    frameon=True,
    loc="center left",
    bbox_to_anchor=(0.80, 0.42),
    fontsize=16,
    borderpad=1.6,
    labelspacing=1.4,
    handletextpad=0.8
)

for legend in [legend1, legend2]:
    for text in legend.get_texts():
        text.set_color("black")
        text.set_fontstyle("normal")
        text.set_fontweight("normal")
    legend.get_frame().set_edgecolor("black")
    legend.get_frame().set_linewidth(1.2)
    legend.get_frame().set_alpha(1.0)

# -----------------------------
# Save
# -----------------------------
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("Wrote:", out_png)
print("Wrote:", out_pdf)