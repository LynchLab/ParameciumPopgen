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
# Use your local paths:
ldhat_path = "poster images/paper figures/data/recombination/pairwise_recombination_table_ldhat.csv"
phi_path   = "poster images/paper figures/data/recombination/phi_full_summary.csv"

# If running in this chat environment with uploaded files, use:
# ldhat_path = "/mnt/data/pairwise_recombination_table_ldhat.csv"
# phi_path   = "/mnt/data/phi_full_summary.csv"

ld = pd.read_csv(ldhat_path)
phi = pd.read_csv(phi_path)

# -----------------------------
# Species labels + colors (your palette)
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
    "tribursaria": "P. tribursaria"
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
    "primaurelia": "#6fa8dc"
}

def norm_species(s):
    s = str(s).strip()
    s = s.replace("Paramecium", "").strip()
    s = s.replace("P.", "").strip()
    s = s.replace(" ", "").replace(".", "")
    return s.lower()

# -----------------------------
# Robust column finder (handles Hudsonâ€™s encoding issue)
# -----------------------------
def find_col(df, must_contain, optional_contain=None):
    must = [x.lower() for x in (must_contain if isinstance(must_contain, (list, tuple)) else [must_contain])]
    opt  = [x.lower() for x in (optional_contain if optional_contain else [])]
    for c in df.columns:
        cl = c.lower()
        if all(m in cl for m in must) and all(o in cl for o in opt):
            return c
    return None

# -----------------------------
# Gene ordering: cox → cob → nadh → rpS → rpL → orf → ymf
# -----------------------------
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

    # anything unexpected goes last, alphabetically
    return (7, g)

# -----------------------------
# Build species order (consistent across panels)
# -----------------------------
phi_species = set(phi["Species"].apply(norm_species)) if "Species" in phi.columns else set()
ld_species  = set(ld["species"].apply(norm_species)) if "species" in ld.columns else set()

species_order = [s for s in [
    "biaurelia","bursaria","caudatum","multimicronucleatum","fokini",
    "quindecaurelia","primaurelia","sexaurelia","tetraurelia","triaurelia","tribursaria"
] if (s in phi_species) or (s in ld_species)]

# ============================================================
# PANEL A: PHI permutation heatmap (gene × species)
# ============================================================
phi2 = phi.copy()
phi2["species_key"] = phi2["Species"].apply(norm_species)
phi2["gene"] = phi2["Gene"].astype(str)

p_col = "PHI_perm"
piv = phi2.pivot_table(index="gene", columns="species_key", values=p_col, aggfunc="min")
piv = piv.reindex(columns=species_order)

sig_thr = 0.05
# 0 = no data, 1 = absent, 2 = present
cat = pd.DataFrame(0, index=piv.index, columns=piv.columns, dtype=int)
cat[(piv >= sig_thr) & (~piv.isna())] = 1
cat[(piv <  sig_thr) & (~piv.isna())] = 2

# keep heatmap readable (main figure)
has_any_signal = (cat == 2).sum(axis=1) > 0
has_enough_data = (cat != 0).sum(axis=1) >= max(3, int(0.25 * len(species_order)))
cat = cat.loc[has_any_signal | has_enough_data]

# ORDER GENES IN YOUR REQUESTED BIOLOGICAL ORDER
cat = cat.loc[sorted(cat.index, key=gene_group_rank)]

# heatmap colors: grey/orange/green
heat_cmap = ListedColormap(["#d9d9d9", "#f4a742", "#2ca25f"])
heat_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], heat_cmap.N)

# ============================================================
# PANEL B: LDhat 4Ner (Hudson vs Wakeley) at smallest θ
# ============================================================
ld2 = ld.copy()
ld2["species_key"] = ld2["species"].apply(norm_species)

# keep gene conversion model
ld2 = ld2[ld2["Model"].str.contains("gene", case=False, na=False)].copy()

theta_col = find_col(ld2, "θ") or find_col(ld2, "theta") or "θ"

# find the two 4Ner columns robustly
hudson_col = find_col(ld2, ["4ner"], ["hudson"]) or find_col(ld2, ["4ner"], ["finite"])
wakeley_col = find_col(ld2, ["4ner"], ["wakely"]) or find_col(ld2, ["4ner"], ["1997"])

if hudson_col is None or wakeley_col is None:
    raise KeyError(
        "Could not identify both 4Ner columns. "
        "Check ld.columns and adjust hudson_col/wakeley_col selection.\n"
        f"Columns seen: {list(ld2.columns)}"
    )

# drop missing and restrict species
ld2 = ld2[ld2[theta_col].notna() & ld2[hudson_col].notna() & ld2[wakeley_col].notna()].copy()
ld2 = ld2[ld2["species_key"].isin(species_order)].copy()

# smallest theta per species
ld_min = (
    ld2.sort_values(["species_key", theta_col])
       .groupby("species_key", as_index=False)
       .first()
)
ld_min = ld_min.set_index("species_key").reindex(species_order).reset_index()

ld_min["label"] = ld_min["species_key"].map(species_labels).fillna(ld_min.get("species", ld_min["species_key"]))
ld_min["color"] = ld_min["species_key"].map(color_map).fillna("#999999")

y_h = ld_min[hudson_col].astype(float).values
y_w = ld_min[wakeley_col].astype(float).values

# ============================================================
# Plot: Heatmap first (A), LDhat second (B)
# ============================================================
fig = plt.figure(figsize=(18.5, 6.8))          # wider canvas
gs = fig.add_gridspec(1, 2, width_ratios=[2.3, 1.0], wspace=0.28)


# --- Panel A: heatmap ---
axA = fig.add_subplot(gs[0, 0])
axA.imshow(cat.values, aspect="auto", cmap=heat_cmap, norm=heat_norm, interpolation="nearest")

axA.set_xticks(range(len(cat.columns)))
axA.set_xticklabels([species_labels.get(sp, sp) for sp in cat.columns], rotation=45, ha="right")
axA.set_yticks(range(len(cat.index)))
axA.set_yticklabels(cat.index)

axA.set_xlabel("Species", fontsize=12)
axA.set_ylabel("Gene", fontsize=12)
axA.set_title("PHI (Pairwise Homoplasy Index) permutation test ", fontsize=13)
axA.text(-0.08, 1.02, "A", transform=axA.transAxes, fontsize=16, fontweight="bold", va="bottom")

# ---- BLACK borders between heatmap cells ----
axA.set_xticks(np.arange(-0.5, cat.shape[1], 1), minor=True)
axA.set_yticks(np.arange(-0.5, cat.shape[0], 1), minor=True)
axA.grid(which="minor", color="black", linestyle="-", linewidth=0.4)
axA.tick_params(which="minor", bottom=False, left=False)

# --- Panel B: LDhat two-method points (NO connecting lines) ---
axB = fig.add_subplot(gs[0, 1])

x = np.arange(len(ld_min))
dx = 0.14
x_h = x - dx
x_w = x + dx

# NO lines between Hudson/Wakeley (requested)

# points: DIFFERENT SHAPES for methods; colors for species
axB.scatter(x_h, y_h, s=85, c=ld_min["color"], marker="o",
            edgecolor="black", linewidth=0.6, zorder=2)
axB.scatter(x_w, y_w, s=85, c=ld_min["color"], marker="s",
            edgecolor="black", linewidth=0.6, zorder=2)

axB.axhline(0, color="black", linestyle="--", linewidth=0.8)

axB.set_xticks(x)
axB.set_xticklabels(ld_min["label"], rotation=45, ha="right")
axB.set_ylabel(r"Population-scaled recombination rate (4$N_e$r)", fontsize=12)
axB.set_title("LDhat gene-conversion estimates", fontsize=13,pad=14)

# ---- Make 0 visible (start slightly below 0) ----
ymax = np.nanmax(np.concatenate([y_h, y_w]))
if not np.isfinite(ymax):
    ymax = 1.0
axB.set_ylim(-0.08, ymax * 1.10)

# integer-like ticks (no 0.5)
axB.set_yticks([0, 1, 2, 3])

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.text(-0.18, 1.02, "B", transform=axB.transAxes, fontsize=16, fontweight="bold", va="bottom")

# -----------------------------
# Legends OUTSIDE (right side)
# -----------------------------
heat_legend = [
    Patch(facecolor="#2ca25f", edgecolor="black", label="Present"),
    Patch(facecolor="#f4a742", edgecolor="black", label="Absent"),
    Patch(facecolor="#d9d9d9", edgecolor="black", label="No data"),
]

method_legend = [
    Line2D([0], [0], marker="o", color="black", markerfacecolor="white",
           markersize=8, linestyle="None", label="Hudson estimator"),
    Line2D([0], [0], marker="s", color="black", markerfacecolor="white",
           markersize=8, linestyle="None", label="Wakeley estimator"),
]

# Make room on the right for legends
plt.subplots_adjust(right=0.80)

fig.legend(handles=heat_legend, title="PHI permutation", frameon=True,
           loc="upper left", bbox_to_anchor=(0.82, 0.86), fontsize=10, title_fontsize=10)

fig.legend(handles=method_legend, title="LDhat 4Ner method", frameon=True,
           loc="upper left", bbox_to_anchor=(0.82, 0.62), fontsize=10, title_fontsize=10)

# Main figure title (PNAS-style)
fig.suptitle(
    "Fig. X. Little evidence for mitochondrial recombination across Paramecium species",
    fontsize=14, y=0.98
)

# -----------------------------
# Save
# -----------------------------
plt.savefig("poster images/paper figures/Images/FigX_PHIperm_LDhat_2methods_main.png",
            dpi=300, bbox_inches="tight")
plt.savefig("poster images/paper figures/Images/FigX_PHIperm_LDhat_2methods_main.pdf",
            bbox_inches="tight")

plt.show()
